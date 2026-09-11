# simphy-x — plan and status

Living document. Newest entries at the bottom of each phase.

## Phase 0 — oracle (done)
- [x] Clone upstream v1.0.2, build with local GSL/MPFR/GMP/SQLite → `upstream/SimPhy/bin/simphy_ref`
- [x] Byte-identical to shipped `simphy_lnx64` on both pipeline command styles
- [x] Profile (gprof) and instrument the bounded-coalescent regime → `DOCS/01-reference-analysis.md`

## Phase 1 — analysis and design (done)
- [x] Full read of main.c / trees.c / num_methods.c / sampling.c / sql_managing.c
- [x] RNG consumption order, numerics, output formats documented
- [x] `DOCS/02-design.md`

## Phase 2 — compat core (byte-exact, single thread) (done 2026-09-11)
- [x] sx_common, sx_rng (mt19937 via GSL), sx_sampling (parser, sampler, stable sort)
- [x] sx_msort (glibc-compatible mergesort)
- [x] sx_tree (SoA trees), sx_stree (BDSA/SSA + outgroup, collapse)
- [x] sx_ltree (SimBDLHTree incl. HGT/GC), CalcProbs, count sampling
- [x] sx_coal (Kahan log series, brent, MPFR chain) ; sx_gtree (MSC, MLC, lineage pool)
- [x] sx_write (all files), sx_db, params/command
- [x] sx_cli (all SimPhy flags), main
- [x] tests/verify_compat.sh green on the 17-case matrix (trees, maps, daughters, stats byte-identical; .params identical
      except the prefix line; .db identical by content; exit codes identical)
- [x] Measured: MSC 10k taxa × 100 loci: RSS 1.58 GB → 10 MB, 2.2 s → 1.4 s; GDL 200 × 20: identical output, 23.6 MB → 4.5 MB
- Not yet supported (Phase 5): -S/-SR/-L/-LR fixed trees, -I config files

## Phase 3 — compat parallel (still byte-exact) (done 2026-09-11, one item deferred)
- [x] record/replay split of the MLC gene-tree simulation (sx_mlc_trace), deferred MPFR chain jobs, pthread pool (sx_tpool)
- [x] per-locus finalisation/writing on workers; DB rows and stats committed in locus order by the main thread
- [x] verify matrix identical for 1/2/4/16 threads (17/17)
- [x] GDL 200 taxa: 20 loci 14.4 s → 3.5 s (16 threads); 100 loci 78 s → 18.8 s, output identical to 1 thread
- [x] MPFR chains distributed per chain (not per locus) over the workers: GDL 200×100 18.8 s → 8.6 s on 16 threads
      (9× vs SimPhy); the sequential skeleton is now 2.5–3 s of that. A parallel probability-row pass was implemented
      but is disabled: the workers are saturated by the chains and the main thread only waited for them.

## Phase 4 — fast mode (done 2026-09-11)
- [x] per-stream RNG (xoshiro256** as gsl_rng_type), stream derivation splitmix64(seed, replicate, locus, tag)
- [x] sx_unif: log-scaled uniformisation (forward transport with adaptive rate/state space, column, endpoint-conditioned
      bridge sampler), Tavaré double series with cancellation check + adaptive-precision MPFR fallback
- [x] unit test tests/unif_test.c: 279 (i,o,t/N) probabilities within 2.6e-10 of 2048-bit MPFR; forward transport
      within 1e-14; bridge quantiles match the conditional CDF (200k samples)
- [x] fast MLC gene-tree path, fast lineage-count probabilities/sampling (sx_ltree.c), fast run loop (sx_run_fast.c)
      with per-job settings snapshot + species-tree clone; ordered commits; thread-count-independent output
- [x] GDL 200×20: 0.96 s single-threaded (SimPhy: 14.4 s); GDL 200×100: 2.45 s / 0.51 s (1 / 16 threads; SimPhy 78 s);
      GDL 500×20 with ≈1000-leaf gene trees: 2.24 s / 0.94 s (SimPhy ≈ 140 s); MSC 10k×100: 0.25 s (16 threads)
- [x] fixed Newick species trees (-S) in both modes; byte-exact vs reference (3 matrix cases)
- [x] statistical equivalence tests vs reference (tests/verify_fast.py): with 200 species trees × 5 loci (20 taxa, GDL
      dup 2.05e-10, Ne 2.9e8) all 11 KS/chi-square comparisons pass (locus leaves/dups/losses, extra lineages, gene
      heights in cu and bl, total lengths, cherries). A 25-tree × 40-loci run showed spurious "failures" whose sign
      flipped between runs: species-tree height varies ±20% (ln:21.25,0.2), so loci are not independent samples and the
      test needs many species trees or a fixed one (-S, Phase 5). Paired diagnostic (identical species and locus trees):
      fast mode with exact vs SimPhy's own count probabilities gives identical gene-tree statistics (all p = 1.0), so the
      lineage-count sampling is equivalent. Pending: fixed-tree paired test; bridge-vs-inversion time-sampler test.
- [x] benchmarks (see 04-verification.md)
- [x] divide-and-conquer bridge for very large Lambda*T: no inversion fallbacks left on the benchmark cases
- [x] exact-CDF tests of the time samplers (tests/bridge_exact_test.c): bridge correct for every coalescence index;
      SimPhy's double-precision attempt produces wrong samples in ~8 % of bounded branches in the cancellation regime
      (documented in 01-reference-analysis.md §7); fixed 100-taxon high-duplication paired test 11/11
- [x] fast mode smoke-tested on all 20 matrix command lines (same exit codes as the reference)

## Phase 4 — fast mode
- [ ] per-stream RNG (xoshiro256** as gsl_rng_type), stream derivation
- [ ] sx_unif: Poisson weights, backward vectors, o_probs, bridge sampler, Tavaré fallback
- [ ] fast MLC gene-tree path; fast count sampling
- [ ] parallel replicate×locus scheduler with ordered output
- [ ] numerics unit tests vs MPFR; statistical equivalence tests vs reference
- [ ] benchmarks 200/500/1000/10000 taxa

## Phase 5 — packaging (done 2026-09-11 except the items marked open)
- [x] Makefile (static deps, `make`, `make tests`, `make check`, `make install`), README, DOCS/03-usage.md, DOCS/04-verification.md
- [x] Fixed Newick species trees (-S) with per-node annotations; configuration files (-I) with the .conf copy
- [ ] Nexus tree files (-SR, -LR) and fixed locus trees (-L): not implemented (settings error), see 03-usage.md
- [x] Pipeline integration: duploss-sim.sh --upto true with --simphy-bin bin/simphy-x gives identical output

## Phase 6 — fidelity and packaging follow-ups (done 2026-09-12)
- [x] SimPhy's stdout text reproduced character for character for verbosity 0–4 (own banner first); per-locus text
      buffered and flushed in locus order in parallel runs; verify_compat.sh checks stdout with SX_CHECK_STDOUT=1 (22/22)
- [x] `make static`: fully static portable binary, byte-identical on the matrix
- [x] GDL memory at scale: 2000 taxa SimPhy 495 MB → compat 33 MB (both abort with SimPhy's MPFR precision error), fast 7 MB
- [ ] still not implemented: -SR/-LR Nexus files, -L fixed locus trees; verbosity 5–6 debugging prints
