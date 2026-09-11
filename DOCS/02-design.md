# simphy-x — design

## 1. Goals

1. **Exactness.** Same model, same distributions, same parameters and the same
   command line as SimPhy 1.0.2. In `compat` mode the output files are
   byte-identical to the reference binary for the same command and seed.
2. **Determinism.** Same command + same seed ⇒ same output, independent of the
   number of threads, in every mode.
3. **Speed and memory.** Multi-core parallelism; O(n) memory instead of O(n²);
   the arbitrary-precision bottleneck removed from the fast path.
4. **Drop-in.** Same flags, same output layout, so `run_simulator.sh`,
   `duploss-sim.sh` and the `concat_*.py`/`reorganize_*.py` scripts work unchanged.

## 2. Two modes

| | `compat` (default) | `fast` |
|---|---|---|
| RNG | one mt19937 stream, SimPhy's draw order | one independent stream per (replicate) and per (replicate, locus), derived from the seed |
| bounded coalescent | SimPhy's double Brent attempt + 512-bit MPFR fallback, op-for-op | uniformisation + endpoint-conditioned path sampling (double precision) |
| lineage-count probabilities | Kahan log-scale alternating series, O(n³) per node | uniformisation, O(M·n) per node |
| output vs SimPhy 1.0.2 | **byte-identical** | same distribution; different sample path |
| parallelism | deferred MPFR chains across cores; per-locus finalisation and writing in parallel; the RNG-consuming skeleton stays sequential | replicates × loci fully parallel |
| expected speed-up | bounded by the sequential skeleton (≈ 5–8× on 16 cores for GDL) | 100× and more on GDL; memory 100–1000× lower at 10k taxa |

Why not one mode: SimPhy's single sequential RNG stream makes byte-identical output
and full parallelism mutually exclusive in general. The parts that *can* be deferred
without changing any random draw are exactly the expensive MPFR chains (their results
never feed back into the number or the order of draws — see 01-reference-analysis §3.4),
so compat mode parallelises those and nothing else.

## 3. Architecture

```
src/
  sx_common.{h,c}   error codes, logging, limits (MAX_LEAVES … from env like SimPhy)
  sx_rng.{h,c}      gsl_rng wrapper: mt19937 (compat) or xoshiro256** custom gsl_rng_type
                    (fast); stream derivation splitmix64(seed, replicate, locus, tag)
  sx_sampling.{h,c} sampling units ("f:", "u:", "ln:", …), parser/printer/sampler,
                    dependency-depth stable sort — port of sampling.c
  sx_msort.{h,c}    glibc-compatible top-down stable mergesort (replaces qsort)
  sx_tree.{h,c}     array-based trees: struct-of-arrays, int32 indices, children via
                    contiguous child lists; species (s), locus (l), gene (g) trees
  sx_stree.{h,c}    NewBDSTree port (BDSA/SSA, outgroup), pre/post-order collapse
  sx_ltree.{h,c}    SimBDLHTree port (dup/loss/HGT/GC, periods, receptors), collapse
  sx_coal.{h,c}     compat numerics: Kahan log-scale series, brent_root, MPFR sampler
  sx_unif.{h,c}     fast numerics: uniformisation vectors, Poisson weights (Fox–Glynn
                    style scaling), endpoint-conditioned path sampler, Tavaré fallback
  sx_gtree.{h,c}    MSC and MLC gene-tree simulation on the lineage-pool stack;
                    CalcProbs and count sampling (compat and fast variants)
  sx_write.{h,c}    Newick writers (%.8lf), mapping/daughters files, params, command
  sx_db.{h,c}       SQLite output (same schema and INSERT text as SimPhy)
  sx_cli.{h,c}      GetSettings port + new options
  sx_run.{h,c}      orchestration: replicate loop, compat skeleton + deferred jobs,
                    fast parallel loop
  main.c
```

### 3.1 Tree representation

All three tree kinds are structs of arrays indexed by node id. Children are stored in
a flat `child[]` array with per-node `(child_off, n_child)`; parent in `anc[]`. Node
payload (times, lengths, kind, paralog, Ne, mu_mult, …) in parallel arrays. Node ids
follow SimPhy's array order (pre-order for the species tree after `CollapseSTree(0)`,
post-order for locus and collapsed gene trees) so that every index printed in an output
file is the same number.

### 3.2 Lineage pool (replaces `l_node::g_nodes`)

One `int32 pool[n_gleaves]` per gene-tree simulation. Processing locus nodes in
post-order, each node owns the segment formed by the segments its children left on
the stack, in child order; coalescences swap-with-last inside the segment exactly as
SimPhy does; the leftover prefix is what the parent sees. Memory O(n_gleaves).

### 3.3 Compat skeleton and deferred MPFR chains

`sx_run` executes SimPhy's loop. Whenever `SimMLCGTree` reaches a bounded locus node
whose double-precision attempt failed, it records a *chain job*:
(n_coals, n_ilin, n_olin, n_gen, gen_length, Ne, epsilon, u[0..n_coals-1]) and stores
placeholder times. Everything downstream that does not depend on the times (pair
choices, gamma multipliers order, topology) proceeds and consumes the RNG exactly as
SimPhy would. As soon as a locus job is produced, its chains are submitted as
individual tasks to the thread pool (one task per bounded node), and the last
finishing chain of a locus triggers the locus finalisation task; this keeps all workers
busy even when a few loci hold most of the MPFR work. A deterministic
finalisation pass then replays the time/branch-length arithmetic per locus in SimPhy's
order (so floating-point sums are bit-identical), and writes the locus outputs. The
locus-level files are independent; the shared files (`l_trees.trees`,
`bounded_locus_subtrees.out`, DB rows) are written from per-locus buffers in locus
order.

`CalcProbsNLineagesLTree` is a post-order tree reduction with per-node deterministic
arithmetic; it is evaluated with a parallel post-order schedule in compat mode as well
(the per-node arithmetic is unchanged, so results are bit-identical).

### 3.4 Fast mode

* Streams: `seed_rep = mix(seed, rep)`, `seed_locus = mix(seed, rep, locus)` with a
  splitmix64 finaliser; each stream is a xoshiro256** exposed through `gsl_rng_type`
  so the GSL distribution samplers (gamma, lognormal, …) are reused unchanged.
* Species-level sampling + species tree on the replicate stream; each locus (params,
  locus tree, all its gene trees, rate heterogeneity, outputs) on its own stream.
  Loci of all replicates are scheduled dynamically over the thread pool; output order
  is fixed by buffering per locus.
* Bounded coalescent per branch (i → o in T): uniformise with Λ = λ_i; backward
  vectors `b_n = Q^n e_o` for n = 0..M; `P(m events | i→o) ∝ Poisson(m;ΛT)·b_m[i]`;
  sample m, sample the bridge forward, event times = sorted uniforms; keep only the
  real (state-changing) events as coalescence times. Branch probability vectors
  `o_probs = Σ_m Poisson(m)·(i_probs Q^m)`. When ΛT exceeds a threshold (default
  4000) the Tavaré series in double with a cancellation check is used; if the check
  fails, adaptive-precision MPFR (rare, kept for completeness).
* Species tree simulation, locus tree simulation, MSC, rate heterogeneity and all
  writers are shared with compat mode.

### 3.5 Command line

All SimPhy 1.0.2 options are accepted with identical semantics. New options
(letters SimPhy does not use): `-XM compat|fast` (`--mode`), `-XT n` (`--threads`,
default: all cores), `-XV` (`--version`). Environment variables `SIMPHY_MAXLEAVES`,
`SIMPHY_MAXIT`, … keep their meaning.

## 4. Verification

1. **Byte-exact compat matrix** (`tests/verify_compat.sh`): a grid of command lines
   (both pipelines' parameterisations, plus outgroup, individuals per species,
   distributions for every sampling unit, HGT/GC with both `-lk`, `-ol`, `-ot`,
   `-om`, `-on`, `-os`, fixed Newick species tree) run through the reference build and
   through simphy-x with 1, 2 and 16 threads; `diff -r` of all files except `.db`
   (compared via `sqlite3 .dump`) and the `.command` name.
2. **Fast-mode statistical equivalence** (`tests/verify_fast.py`): many replicates
   of small/medium scenarios in both tools; compare distributions of gene-tree
   height, total length, number of extra lineages, number of locus-tree leaves,
   duplication/loss counts, coalescence-time quantiles and topology summaries
   (Kolmogorov–Smirnov / chi-square with Bonferroni).
3. **Unit tests for the numerics**: uniformisation vs MPFR-512 Tavaré on a grid of
   (i, o, T/Ne): relative error < 1e-10 for both P_{i,o} and the CDF of the first
   waiting time (Monte-Carlo check of the bridge sampler against the MPFR CDF).
4. **Benchmarks** (`bench/`): time and RSS for 200/500/1000/10000 taxa, MSC and GDL,
   1 vs 16 threads, versus the reference.
