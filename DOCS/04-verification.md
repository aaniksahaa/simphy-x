# simphy-x — verification

All numbers below were measured on the development machine (16 cores, glibc 2.39,
GSL 2.7.1, MPFR 4.2.1) on 2026-09-11. Scripts live in `tests/`.

## 1. Oracle

`upstream/SimPhy/bin/simphy_ref` is SimPhy v1.0.2 built from the upstream tag with the
system libraries. Its output is byte-identical to the official SimPhy 1.0.2 Linux release `simphy_lnx64` on both
command styles we use (multispecies coalescent and duplication/loss pipelines).
Build it with `make oracle` (after `git submodule update --init`).

## 2. Byte-exact compat mode — `tests/verify_compat.sh`

Runs 28 command lines through the oracle and through `simphy-x` (`-XM compat`, default) with
the requested thread counts, and compares:

* every file under the output directory except `*.command`, `*.params`, `*.db`: `diff -r`
  (species/locus/gene trees, `.mapsl`/`.maplg` mappings, `bounded_locus_subtrees.out`,
  `stats.txt`);
* `*.params`: identical except the "Output files prefix" line (it contains the directory name);
* `*.db`: identical table contents (SQLite files themselves embed the library version);
* exit codes (error paths included);
* with `SX_CHECK_STDOUT=1`, the program text after the banner separator (verbosity 0–4;
  the output-directory name is normalised).

Cases: MSC pipelines (fixed and sampled parameters, 500 taxa × 50 loci), outgroup + several
individuals per species + time units + node labels + stats, SSA species trees (with/without
outgroup), sampled tree origin, all rate-heterogeneity layers, GDL pipelines (mappings, labels,
two individuals), HGT with both receptor-choice modes, gene conversion, genome-wide priors with
dependent sampling units, minimum-leaves settings (error path), generation-time units, fixed
Newick species trees (plain, with per-node Ne and individual counts, with GDL+HGT), Nexus
species-tree files (`-SR`, plain and with `[&pop_size, n_ind, u_mult]` comments, verbosity 4),
fixed locus trees (`-L`, duplication and loss nodes, per-leaf individual counts, gene-tree
rate heterogeneity, labels), Nexus locus-tree files (`-LR`, two replicates over three trees,
generation time 2). The inputs live in `tests/inputs/`.

Result: **28/28 identical (files and stdout) with 1, 4 and 16 threads**, also with the
dynamically linked build `bin/simphy-x-dynamic`. Larger checks: GDL 200 taxa × 20 and
× 100 loci identical to the oracle (1 and 16 threads), MSC 10 000 taxa × 100 loci identical.

## 3. Fast-mode numerics — `tests/unif_test.c`

* 279 (i, o, t/Ne) combinations with Ne = 2.9e8, i ≤ 120, o ≤ 60, t/Ne from 1e-6 to 5:
  uniformisation vs 2048-bit MPFR Tavaré series, worst relative error 2.6e-10.
* forward transport of lineage-count distributions (fmax = 40): worst relative error 1e-14.
* endpoint-conditioned bridge sampler: 200 000 paths for (i=30, o=12, T/Ne=1e-3); the
  empirical quantiles of the first waiting time agree with the exact conditional CDF within
  Monte-Carlo tolerance.

`tests/bridge_vs_inversion.c` compares the bridge sampler with SimPhy's own sequential CDF
inversion (double attempt + 512-bit MPFR redo, the exact compat logic) on fixed parameters
using Kolmogorov–Smirnov statistics on every waiting time and on the total (5000 paths each):
4 of 5 configurations agree (max KS ≤ 0.038, critical 0.039); (i=40, o=20, T/Ne=1e-3) does
not (KS 0.073). `tests/bridge_exact_test.c` then tests each sampler separately against the
**exact** marginal CDF of the j-th coalescence time P(τ_j ≤ x) = Σ_{k≤i−j} P_{i,k}(x)P_{k,o}(T−x)
/ P_{i,o}(T) at 2048 bits: the bridge sampler passes for every j (KS ≤ 0.016, n = 5000–20000),
SimPhy's sampler fails for every j; splitting SimPhy's paths by the code path that produced
them shows the MPFR redo is correct (KS 0.040) and the accepted double-precision attempts are
wrong (KS 0.98, ~8 % of paths). `tests/simphy_cdf_check.c` shows SimPhy's MPFR CDF formula is
exact and its double CDF is ±inf in that regime. See 01-reference-analysis.md §7.

## 4. Fast mode vs SimPhy statistically — `tests/verify_fast.py`

Same command line in both modes, many loci; per-locus (leaves, duplications, losses) and
per-gene-tree (extra lineages, height in coalescent units, total length, cherries, Newick
height/length in substitution units) statistics compared with two-sample KS and chi-square
tests, Bonferroni-corrected at 1 %.

Important design point: statistics that depend on replicate-level random effects (species
tree height and topology, lineage-specific rate multipliers drawn once per replicate from a
gamma with a lognormally sampled shape) make the pooled loci non-independent, so the test must
either sample many species trees or fix the species tree (`--stree`) and the replicate-level
rates (`--simple-rates`).

Results (GDL, dup 2.05e-10, loss 0.1×, Ne 2.9e8, 1000 loci per side):
* 200 species trees × 5 loci: 11/11 pass.
* fixed 20-taxon tree, no replicate-level rate heterogeneity, 20 × 50 loci: 11/11 pass;
  means agree to < 1 % (e.g. gene total length 14.61 vs 14.62, extra lineages 104.6 vs 102.9).
* MSC only (no duplications), fixed 20-taxon tree, 20 × 50 loci: 8/8 pass (extra lineages 57.33 vs 57.34).
* paired diagnostic with identical species *and* locus trees (fast mode with exact vs with
  SimPhy's own count probabilities): all statistics identical (p = 1.0), i.e. the
  lineage-count sampling is equivalent to SimPhy's.
* SimPhy's double-precision count probabilities were logged on a run (24 308 evaluations) and
  compared with exact values: the ones used for sampling are exact; in the bottom-up pass
  1 433 evaluations return 0 and 1 866 have > 1e-6 relative error, all at probabilities below
  1e-11, i.e. without practical effect.

## 5. Performance (wall time / peak RSS)

| case | SimPhy 1.0.2 | simphy-x compat, 16 thr | simphy-x fast, 1 thr | simphy-x fast, 16 thr |
|---|---|---|---|---|
| GDL 200 taxa × 20 loci | 14.6 s / 24 MB | 3.5 s / 12 MB | 0.96 s / 29 MB | — |
| GDL 200 taxa × 100 loci | 78 s | 8.6 s / 20 MB | 1.97 s / 19 MB | 0.42 s / 222 MB |
| GDL 500 taxa × 20 loci, dup ≈ 3 copies (≈ 1000-leaf gene trees) | ≈ 140 s (extrapolated) | 13.5 s / 46 MB | 1.88 s / 34 MB | 0.76 s / 158 MB |
| MSC 10 000 taxa × 100 loci | 2.2 s / **1.58 GB** | 0.78 s / 23 MB | — | 0.25 s / 117 MB |
| GDL 2000 taxa × 3 loci (dup 6e-11) | aborts after 5.2 s / 495 MB ("MPFR out of precision") | aborts identically, 33 MB; with `-XR 1`: completes in 9.9 s / 44 MB (2 inversions needed > 512 bits) | — | completes, 7 MB |

Compat mode's remaining sequential part (species/locus trees, count sampling, the
double-precision attempts, pair choices) is 2.5–3 s of the 200×100 run; everything else is the
MPFR chains, which are distributed per chain over the workers. Fast mode has no sequential
part beyond the species trees; its 16-thread numbers on 20–100 loci are limited by the
per-locus granularity (few loci per thread).

## 6. Pipeline integration

Our duplication/loss pipeline script (`duploss-sim.sh --upto true --simphy-bin bin/simphy-x`) produces a
`true-genetrees/` and `species-tree/` layout identical to the run with the shipped binary
(checked with `diff -r` on a 10-taxa × 5-loci × 2-replicate run).
