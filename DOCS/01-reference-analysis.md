# SimPhy 1.0.2 — reference analysis

This document records what the reference simulator does, in enough detail to
reproduce its output byte-for-byte, and where its time and memory go.

## 1. Provenance and oracle

* Shipped binary: `ref/*/simphy_lnx64` (md5 `3969a8370c84fc379d62f405ec8f7867`,
  GCC 4.8.2, March 2016, GSL/MPFR/GMP/SQLite statically linked, libc/libm dynamic).
* Source: upstream tag `v1.0.2` of https://github.com/adamallo/SimPhy, cloned to
  `upstream/SimPhy`. Built here with the system GSL 2.7.1 and locally extracted
  MPFR 4.2.1 / GMP 6.3.0 / SQLite 3.45.1 (`deps/`), no root needed.
* **The fresh build reproduces the shipped binary byte-for-byte** on both pipeline
  argument styles (normal MSC and GDL). Hence `upstream/SimPhy/bin/simphy_ref` is a
  valid, compilable, instrumentable oracle.

## 2. Measured behaviour (this machine, 16 cores)

| case | wall | max RSS |
|---|---|---|
| normal, 1000 taxa × 100 loci | 0.20 s | 20 MB |
| normal, 10000 taxa × 100 loci | 2.3 s | **1.5 GB** |
| GDL, 200 taxa × 50 loci (dup 2.0e-10, loss 0.1×) | **50.6 s** | 41 MB |

* The GDL (duplication/loss) path is the time bottleneck: ~1 s per locus at 200 taxa,
  which extrapolates to hours for the SIM200/SIM500 datasets and is unusable at 1000+.
* The normal path is the memory bottleneck: every locus-tree node allocates an array of
  `n_gleaves` gene-node pointers (`NewLNodes`), i.e. O(nodes × gene leaves) = 20001 ×
  10001 × 8 B ≈ 1.6 GB at 10k taxa. It is quadratic.

### Profile of the GDL run (gprof, 20 loci)

~88 % of the time is GMP/MPFR arithmetic (`__gmpn_*`, `mpfr_exp_2`, `mpfr_mul`, …)
called from `MPFRSampleCDFCoalTimeMLCFromXtoYLineages` (115 857 evaluations inside
Brent root-finding, 512-bit precision). ~4 % is `KahanLogscaleProbCoalFromXtoYLineages`
(7.6 M `log1p` calls) from `CalcProbsNLineagesLTree`. Everything else is negligible.

### Regime of the bounded-coalescent sampler (instrumented run, 200 taxa)

| path | calls | i_lin | o_lin | T/Ne (median, max) | λ_i·T (median, max) |
|---|---|---|---|---|---|
| double attempt (Kahan log-scale + Brent) | 1309 | ≤ 63 | 1–10 | 1.93, 5.7 | 54, 2.3e3 |
| MPFR fallback (512 bit) | 2296 | 10–63 | 18–28 | **9e-4, 4e-3** | **0.64, 8.2** |

The MPFR fallback is triggered when the branch is very short relative to Ne and many
lineages exit it: all exponentials in Tavaré's alternating series are ≈ 1 and the
binomial-sized coefficients cancel to a tiny result, so doubles fail. In that same
regime the number of expected uniformised events λ_i·T is ≤ 10, which makes the
positive-term uniformisation series (Section 6) essentially free.

## 3. Program flow and random-number consumption order

One `gsl_rng` (`gsl_rng_default` = mt19937, `gsl_rng_set(r, seed)`) is shared by
everything and consumed strictly sequentially. Reproducing SimPhy exactly means
reproducing this order of draws.

```
GetSettings (no RNG)
for replicate s = 1..RS:
  sample_distr(17 species-level sampling units)     # stable-sorted by dependency depth
  NewBDSTree(...)                                    # species tree (BDSA or SSA)
  CollapseSTree(pre-order)
  if alpha_s != 0: Rateheter_lineagespec -> one gsl_ran_gamma per s_node, array order
  write s_tree.trees ; DB species row
  for locus l = 1..RL:
    sample_distr(lb,ld,lt,lgc,lalpha_g)              # pointer-dependent: no draws unless distributions
    if alpha_l != 0: gamma_l = gsl_ran_gamma(alpha_l, 1/alpha_l)
    if any of lb,ld,lt,lgc != 0: SimBDLHTree (locus tree, retries until valid)
         ReindexSTree(post-order); write .mapsl if map
         gene tree allocated with 2*n_gleaves-1 nodes
    else: locus tree = copy of species tree (once)
    write l_trees.trees line; daughters file line; DB locus row
    for gene tree g = 1..RG:
      sample_distr(alpha_g)
      SimMSCGTree (no dup/transfer/gc) or SimMLCGTree (bounded, multilocus coalescent)
      if alpha_g != 0: Rateheter_GTbranchspec -> one gsl_ran_gamma per g_node, array order
      write g_treesNNN.trees ; .maplg if map ; DB gene row
```

### 3.1 Species tree, `NewBDSTree` (BDSA, `-sl` given)

With `-sl N` and `-st T` both given (both pipelines), the "TMRCA given" branch runs:
`n1 = gsl_rng_uniform_int(leaves-1)+1` splits N into two BDSA sub-processes (special
cases: leaves==2 → no loop; n1==1 → single process with N-1 leaves). For each
sub-process: `n-1` uniforms `gsl_rng_uniform_pos` → quantile function of speciation
times (Hartmann/Gernhard), `qsort` ascending (epsilon comparator `Compare_DBL`), then
reconstruct backwards choosing pairs with `gsl_rng_uniform_int(avail)` (second pick
only if `avail>1`). With `-so ratio` an outgroup leaf (sp_index **0**) is attached at
height `(1+ratio)·T`. Leaves get `sp_index = j+1`. `time==0` (only `-sl`) samples the
origin from the conditional CDF instead. Without `-sl` the SSA algorithm runs.

### 3.2 Locus tree, `SimBDLHTree`

Pre-order walk over species nodes 1..n-1; for every locus lineage present in the
parent species branch: `-log(U)/w_prob` waiting time; if it exceeds the branch, one
child l_node; else a birth–death(+HGT+GC) SSA along the branch: `uniform_int(avail)`
(only if avail>1) picks the lineage, `U=gsl_rng_uniform` picks the event
(`>dtgc_prob` duplication, `>tgc_prob` loss, `>gc_prob` transfer donor, else GC donor),
then `-log(U)/(w_prob*avail)`. Duplication child 1 gets `paralog=next_paralog++`,
`fmax_nlin=1, n_olin=1` (the bounded subtree root). Retries whole tree while
`lt_true_leaves<min_lleaves || MAX_LEAVES reached || diffs<min_lsleaves || gleaves<=1`.
Transfers/GCs are completed afterwards with a period sweep (`qsort` of periods by
time, `Compare_periods`) and receptor choice (`ChooseLNodePeriod(gsl_rng_uniform_pos)`
or `uniform_int`). Finally `CollapseResizeLTree(post-order)` and, if any dup/transfer/gc,
`CalcProbsNLineagesLTree` (no RNG).

Safety check before every trial: `exp((b-d)·H)·n_leaves >= MAX_LEAVES` → error
(`SIMPHY_MAXLEAVES` env var raises the cap; defaults MAX_LEAVES=1000, MAX_IT=10000).

### 3.3 Gene tree, `SimMSCGTree` (no events)

Post-order over locus nodes. Node's active lineages = its gene leaves (leaf) or the
concatenation of the children's leftovers in child order. Loop: `t = -log(U)·2Ne/(k(k-1))`;
if the branch end is passed (and not root) stop; else pick `off1 = list[(int)(U·avail)]`,
swap-with-last, `off2 = list[(int)(U·avail)]` (only if avail>1, else 0), replace slot
with the new internal node (index = n_gleaves + running counter). Branch lengths
`bl += (min(child.n_gen, branch top) - anc.n_gen)·mu·mu_mult` and leftovers get
`bl += (…)·mu` when handed to the parent. Extra lineages counted.

### 3.4 Gene tree, `SimMLCGTree` (bounded / multilocus coalescent)

1. `MatchTreesMLC`: gene leaves numbered in locus post-order; `n_olin=1` for bounded
   roots (child 1 of DUP/TRFR/GC with `fmax_nlin>0`), else 0.
2. `SampleNLineagesCDFLTree` (top-down): for a node with `n_olin>0`: `u=gsl_rng_uniform`,
   walk `kt = n_olin..ktmax` accumulating `P(kt→n_olin | t)·i_probs[kt]/o_probs[n_olin]`
   (Kahan) until `>= u`; then split kt between two children with `SampleNDoubles`
   over `i_combprobs` (one `gsl_rng_uniform_pos`), or give it to the single child.
3. Per locus node in post-order: unbounded (`n_olin==0`): exponential waiting times,
   pre-drawn (`-log(U_pos)·2Ne/(k(k-1))`) while inside the branch. Bounded: for
   `j<n_coals`: `SampleCoalTimeMLCFromXtoYLineages(avail, n_olin, remaining, Ne, eps,
   U_pos)` (double Kahan-log series + Brent), break early when the branch is exhausted;
   if the last time `<=0` or the branch is exhausted → **redo all `n_coals` with
   `MPFRSampleCoalTimeMLCFromXtoYLineages(…, U_pos, 512 bits)`** (fresh uniforms).
   Then the pair-picking loop exactly as in MSC (`(int)(U·avail)` twice).
4. `CollapseGTree(post-order)` only when mapping or labels are requested. Node array
   order determines which gamma multiplier each branch receives in
   `Rateheter_GTbranchspec`, so it affects output.

### 3.5 Numerics that must be replicated exactly for byte-identical output

* `epsilon` is a `float` (default 1e-6f) promoted to double.
* `brent_root` (num_methods.c:53): modified Brent (refuses to return the bracket ends).
* `KahanLogscaleProbCoalFromXtoYLineages`, `KahanLogscaleAdd`, `LogscaleSub`,
  `LogscaleFact` (255-entry table, Stirling above), `biomcmc_log1p` (Chebyshev
  series, falls back to `gsl_log1p`), `SKahanSum`, `VKahanSum`.
* MPFR: `MPFROriginalProbcoalFromXtoYLineages` and `MPFRSampleCDF…` at 512 bits,
  `MPFR_RNDN`, result via `mpfr_get_ld` then stored into a `double`. MPFR is correctly
  rounded per operation, so the op sequence fixes the bits independently of version.
* `qsort` with an epsilon comparator: glibc (2.39 here) uses a top-down stable
  mergesort; simphy-x re-implements that exact algorithm (`sx_msort`) rather than
  calling `qsort`, so the behaviour no longer depends on libc.
* All exp/log go through libm (dynamic in the reference too); bit-exactness is
  therefore "same machine/libm" exactness, identical to the reference's own guarantee.

## 4. Output files and formats

Directory `-o NAME/` (chdir'd into), per replicate a zero-padded dir `%.*d`
(`count_intdigits(RS)` digits). In it:

* `s_tree.trees`: `WriteSNodesFileGen`: leaves `%d:%.8lf` (sp_index; names only for
  parsed trees), internal `(...):%.8lf`, root `(...);\n`. `-ot 1` multiplies by
  `gtime_mult·gen_time`; `-ol 1` prints `)%d:` with the node index.
* `l_trees.trees`: one line per locus, leaves `%d_%d:%.8lf` (sp_index_paralog), lost
  leaves `Lost-%d_%d` (species post-order index), `Rtransf-`, `Rgc-`.
* `g_trees%.*d.trees` (digits = `count_intdigits(RL)`): leaves `%d_%d_%d:%.8lf`
  (sp_index_paralog_replica, bl), internal `):%.8lf` (`)%d:` with `-ol`), root `);\n`.
* `bounded_locus_subtrees.out` (`-on 1`): per locus, daughters of DUP/TRFR/GC nodes.
* `%.*d.mapsl`, `%.*dl%.*dg.maplg` (`-om 1`): tab-separated mappings (formats in
  `WriteMappingSL/LG`).
* `NAME.command` (argv joined with spaces, trailing space), `NAME.params`
  (`PrintGlobalSettings`), `NAME.db` (SQLite, 3 tables), `stats.txt` (`-os 1`).

## 5. Memory model (why it is quadratic)

`l_node::g_nodes` is an array of `n_gleaves` pointers per locus-tree node, allocated
in `NewLNodes`, `CollapseLTree`, `RefineLNodes`. Only the concatenation of the
children's leftover lineages is ever stored in it, so the total live content is at
most `n_gleaves` pointers at any time. Because nodes are processed in post-order,
the lists form a stack: a node's list is exactly the contiguous top-k segments pushed
by its children. simphy-x therefore uses one pool of `n_gleaves` indices and per-node
(offset, count) pairs — same order, same swap-with-last semantics, O(n) memory.

## 6. Why MPFR can go: the bounded coalescent as an endpoint-conditioned CTMC

Inside a bounded subtree the lineage count along a branch is a pure-death chain with
rates λ_k = k(k−1)/(2Ne), conditioned on going from i lineages to o lineages in
time T. SimPhy samples successive waiting times by inverting the conditional CDF,
whose closed form (Tavaré 1984) is an alternating series with binomial-sized terms.

Uniformisation (Λ = λ_i, DTMC Q with jump probability λ_k/Λ) gives every needed
quantity as a **positive** series: P_{k,o}(t) = Σ_m Poisson(m; Λt)·(Q^m)[k,o], and
the whole vector o_probs = i_probs·e^{Qt} = Σ_m Poisson(m)·(i_probs Q^m) in O(M·n)
instead of O(n³) alternating evaluations. Endpoint-conditioned paths are sampled
exactly with the Hobolth–Stone (2009) construction: draw the number of uniformised
events m ∝ Poisson(m; ΛT)·(Q^m)[i,o], then the event times are i.i.d. uniform on
[0,T], and the state path is a forward-sampled bridge using the backward vectors
Q^{m−j}e_o. No root finding, no cancellation, no arbitrary precision. The instrumented
regime above (ΛT ≤ 8 where MPFR was needed, ≤ 2300 elsewhere) makes this cheap; for
very large ΛT the well-conditioned Tavaré series is used, with adaptive-precision MPFR
kept only as a last-resort fallback.

## 7. A numerical defect in SimPhy 1.0.2 (found while validating fast mode)

For every bounded locus branch SimPhy first samples the coalescence waiting times with a
double-precision Kahan log-scale CDF and Brent's method, and only redoes the branch with
512-bit MPFR when the last time is ≤ 0 or the branch is exhausted (`SimMLCGTree`,
trees.c:5524, the `if(*(sampled_ngens+n_coals-1)<=0 || …)` test). In the cancellation regime
(short branch relative to Ne, many lineages, e.g. i=40 → o=20 in T/Ne = 1e-3) the double CDF
evaluates to ±inf (`tests/simphy_cdf_check.c`), Brent still returns a value inside (0, T), and
the acceptance test does not catch it. `tests/bridge_exact_test.c` compares SimPhy's samples
with the exact marginal CDF of the j-th coalescence time (2048-bit MPFR): the paths produced by
the MPFR redo are correct (KS 0.040, n = 734), the paths accepted from the double attempt are
essentially arbitrary (KS 0.98, n = 66, i.e. ~8 % of the branches in that regime), and the
pooled sample fails for every j. SimPhy's MPFR CDF formula itself is correct (it matches the
exact CDF to 8 digits in all tested configurations).

Consequences: compat mode reproduces this behaviour exactly (it is part of SimPhy's output);
fast mode samples the exact conditional distribution (validated against the same exact CDFs,
KS ≤ 0.016 with n = 5000–20000 for every j). At the gene-tree level the effect is small: the
paired fast-vs-compat comparisons of extra lineages, heights and lengths pass, and the mean
total coalescence time in the affected branches is unchanged (0.000954 vs 0.000952 Ne).
