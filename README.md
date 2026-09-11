# simphy-x

A from-scratch, parallel re-implementation of the phylogenomic simulator
[SimPhy](https://github.com/adamallo/SimPhy) 1.0.2 (species trees, gene-family/locus trees with
duplication, loss, transfer and gene conversion, gene trees under the multispecies and bounded
multilocus coalescent, rate heterogeneity).

* **compat mode (default)**: byte-identical output to SimPhy 1.0.2 for the same command line
  and seed, with any number of threads (verified on a 20-case matrix); O(n) memory instead of
  O(n²) (10 000 taxa: 1.6 GB → 10 MB); the arbitrary-precision work of duplication/loss
  simulations runs in parallel (200 taxa × 100 loci: 78 s → 8.6 s on 16 cores).
* **fast mode** (`-XM fast`): same model, independent random streams per locus, uniformisation
  instead of 512-bit MPFR root finding; 40× faster than SimPhy single-threaded and 180× on
  16 cores for the duplication/loss pipeline (200 taxa × 100 loci: 78 s → 0.42 s), statistically
  indistinguishable output, and it samples the bounded coalescent exactly where SimPhy 1.0.2's
  double-precision shortcut occasionally does not (DOCS/01-reference-analysis.md §7).

```
make
bin/simphy-x -h
bin/simphy-x <SimPhy options> [-XM compat|fast] [-XT threads]
```

Documentation: `DOCS/00-plan.md` (status), `01-reference-analysis.md` (what SimPhy does and
why it is slow), `02-design.md`, `03-usage.md`, `04-verification.md`.
Tests: `tests/verify_compat.sh` (byte-exact matrix against the reference build),
`tests/verify_fast.py` (statistical equivalence), `tests/unif_test.c`,
`tests/bridge_vs_inversion.c`, `tests/compat_prob_test.c` (numerics).
Reference source: `upstream/SimPhy` (tag v1.0.2); reference binaries in `ref/`.
