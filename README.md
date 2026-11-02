# SimPhy-X

**A scalable, parallel, SimPhy-compatible simulator of species, locus and gene trees.**

[![License: GPL v2+](https://img.shields.io/badge/license-GPL--2.0--or--later-blue.svg)](LICENSE)
[![Based on SimPhy 1.0.2](https://img.shields.io/badge/based%20on-SimPhy%201.0.2-green.svg)](https://github.com/adamallo/SimPhy)

SimPhy-X is a modified and re-engineered version of
[SimPhy 1.0.2](https://github.com/adamallo/SimPhy) by Diego Mallo, Leonardo de Oliveira Martins
and David Posada ([Mallo et al., 2016](https://doi.org/10.1093/sysbio/syv082)). It simulates the
same models, which are:

* species trees (birth–death),
* locus trees (gene duplication, loss, horizontal transfer and gene conversion),
* gene trees (the multispecies and bounded multilocus coalescent), with
* substitution-rate heterogeneity.

It accepts **SimPhy's command line unchanged** and writes the same files. The difference is
scale. SimPhy-X uses every core, needs linear rather than quadratic memory, and can sample
much faster without arbitrary-precision arithmetic.

> SimPhy-X is an independent project. It is not affiliated with or endorsed by the SimPhy
> authors. Credit for the SimPhy model, its design and its original implementation belongs to
> them. Please cite their paper as well as SimPhy-X (see [Citation](#citation)).

## Two modes

| | `compat` (default) | `fast` (`-XM fast`) |
|---|---|---|
| Output vs SimPhy 1.0.2 | **byte-identical** for the same command and seed | same model and distributions; a different random sample |
| Threads | any; the output never depends on the thread count | any; the output never depends on the thread count |
| Memory | O(n) instead of SimPhy's O(n²) | O(n) |
| Bounded coalescent | SimPhy's algorithm, operation for operation (512-bit MPFR) | uniformisation and exact endpoint-conditioned sampling (double precision) |
| Use it for | reproducing published SimPhy datasets | new and large datasets |

## Scalability

Wall time and peak memory on a 16-core Linux machine. All numbers come from
[DOCS/04-verification.md](DOCS/04-verification.md). GDL means gene duplication and loss; MSC
means the multispecies coalescent without duplications.

| Simulation | SimPhy 1.0.2 | SimPhy-X compat, 16 threads | SimPhy-X fast, 1 thread | SimPhy-X fast, 16 threads |
|---|---|---|---|---|
| GDL, 200 taxa × 20 loci | 14.6 s / 24 MB | 3.5 s / 12 MB | 0.96 s / 29 MB | — |
| GDL, 200 taxa × 100 loci | 78 s | 8.6 s / 20 MB (**9×**) | 1.97 s / 19 MB (**40×**) | 0.42 s (**185×**) |
| GDL, 500 taxa × 20 loci (≈1000-leaf gene trees) | ≈140 s (extrapolated) | 13.5 s / 46 MB | 1.88 s / 34 MB | 0.76 s |
| MSC, 10 000 taxa × 100 loci | 2.2 s / **1.58 GB** | 0.78 s / **23 MB** | — | 0.25 s |
| GDL, 2000 taxa × 3 loci | aborts ("MPFR out of precision") | aborts identically; with `-XR 1` completes in 9.9 s / 44 MB | — | completes, 7 MB |

* **Memory.** SimPhy allocates a gene-node array for every locus-tree node, which is quadratic
  in the number of taxa. SimPhy-X keeps trees in flat arrays, so 10 000 taxa need 23 MB
  instead of 1.58 GB, about 70× less.
* **Compat-mode parallelism.** SimPhy's single random stream forces the random draws into one
  sequential pass. SimPhy-X runs that pass on one thread and spreads the expensive MPFR work,
  file writing and statistics across all the others, without changing a single draw.
* **Fast mode.** Every (replicate, locus) pair gets its own random stream, so loci run fully
  in parallel. The bounded coalescent is sampled exactly in double precision, which removes
  the MPFR bottleneck.
* **Correctness.** Fast mode also avoids a numerical defect in SimPhy 1.0.2: about 8% of
  bounded-coalescent branch samples in the cancellation regime are accepted from a
  double-precision computation that returns ±inf. This is documented in
  [DOCS/01-reference-analysis.md §7](DOCS/01-reference-analysis.md). Compat mode reproduces
  the defect on purpose, because it must match SimPhy byte for byte.

## Installation

### Prebuilt binary (x86-64 Linux)

`bin/simphy-x` is a fully static executable. It needs no libraries and gives the same output on
every machine.

```bash
curl -L -o simphy-x https://github.com/aaniksahaa/simphy-x/raw/main/bin/simphy-x
chmod +x simphy-x
./simphy-x -XV
```

### Build from source

You need `gcc`, `make` and the GSL development package (`libgsl-dev` on Debian/Ubuntu). The
MPFR, GMP and SQLite headers and static libraries for x86-64 are bundled in `deps/`.

```bash
git clone https://github.com/aaniksahaa/simphy-x.git
cd simphy-x
make                                  # bin/simphy-x, fully static and portable
make install PREFIX=$HOME/.local      # optional: installs $HOME/.local/bin/simphy-x
```

To link against the system's glibc instead, run `make dynamic`, which builds
`bin/simphy-x-dynamic`. To build against system-wide MPFR, GMP and SQLite, add `DEPS=/`.

## Quick start

Any SimPhy 1.0.2 command works as it is: replace `simphy` with `simphy-x`.

**Multispecies coalescent.** This simulates 10 replicates of a 100-taxon species tree with
1000 gene trees each, and the output is identical to SimPhy's:

```bash
simphy-x -rs 10 -rl f:1000 -rg 1 -sl f:100 -sb f:0.000001 -st f:1000000 \
         -sp f:10000 -su f:0.00001 -sg f:1 -cs 22 -o msc_100
```

**Gene duplication and loss, fast mode on all cores:**

```bash
simphy-x -rs 20 -rl f:1000 -rg 1 -sl f:200 -sb f:0.000001 -st ln:21.25,0.2 -so f:1 \
         -sp f:290000000 -su ln:-21.9,0.1 -hs ln:1.5,1 -hl ln:1.551533,0.6931472 -hg ln:1.5,1 \
         -lb f:2.046e-10 -ld f:2.046e-11 -cs 9644 -o gdl_200 -XM fast
```

**Fewer threads.** Add `-XT 8` to use 8 threads. The output is the same for every thread
count.

**Fixed species tree:**

```bash
simphy-x -s "(((A:1000000,B:1000000):500000,(C:800000,D:800000):700000):300000,E:1800000);" \
         -rs 2 -rl f:100 -rg 1 -sp f:20000 -su f:0.00001 -cs 61 -o fixed
```

Each run writes `<prefix>/<replicate>/s_tree.trees`, `l_trees.trees` and `g_trees*.trees`.
Optional outputs follow the SimPhy flags: mappings (`-om 1`), an SQLite database (`-od 1`),
sampled parameters (`-op 1`), bounded locus subtrees (`-on 1`), per-replicate statistics
(`-os 1`) and others. [DOCS/cookbook.sh](DOCS/cookbook.sh) collects more recipes, including
horizontal transfer, gene conversion, genome-wide priors, Nexus inputs, fixed locus trees and
configuration files.

## Command line

```
simphy-x [SimPhy 1.0.2 options] [-XM compat|fast] [-XT threads] [-XR 0|1]
simphy-x -h          # full option list
simphy-x -XV         # version
```

All SimPhy 1.0.2 options keep their meaning and defaults; see `simphy-x -h` and the
[SimPhy manual](https://github.com/adamallo/SimPhy/wiki). SimPhy-X adds:

| Option | Meaning |
|---|---|
| `-XM compat` | Default. Output byte-identical to SimPhy 1.0.2 for the same command and seed. |
| `-XM fast` | Same model; an independent random stream per (replicate, locus); no arbitrary-precision arithmetic. Deterministic for a given seed, but not SimPhy's sample path. |
| `-XT n` | Number of threads. The default is all available cores. Never changes the output. |
| `-XR 1` | Compat mode only. Where SimPhy 1.0.2 aborts with "MPFR out of precision" (bounded subtrees with hundreds of lineages), continue with higher precision. Output up to that point is identical. |
| `-XV` | Print the version and exit. |

Long forms `--mode`, `--threads`, `--version` and `--help` are also accepted. SimPhy's
environment variables (`SIMPHY_MAXLEAVES`, `SIMPHY_MAXIT`, `SIMPHY_MAXNAME`,
`SIMPHY_FLOATPREC`) keep their meaning. For large duplication rates, raise
`SIMPHY_MAXLEAVES`.

### Which mode should I use?

* Use **compat** to regenerate a dataset that was made with SimPhy 1.0.2, or whenever results
  must match SimPhy byte for byte.
* Use **fast** for new datasets, large taxon counts, and duplication/loss models. Its
  distributions are statistically indistinguishable from SimPhy's (see below), and it
  completes cases that SimPhy cannot run at all.

## Verification

* **Byte-exact compat mode.** `tests/verify_compat.sh` runs 28 command lines through SimPhy
  1.0.2, built from the upstream tag, and through SimPhy-X. The cases cover every input form
  (`-S`, `-SR`, `-L`, `-LR`, `-I`), HGT, gene conversion, outgroups, all rate-heterogeneity
  layers and error paths. It compares every output file and the program text. Result: 28/28
  identical with 1, 4 and 16 threads.
* **Fast mode.** `tests/verify_fast.py` compares locus-tree and gene-tree statistics between
  the two modes with KS and chi-square tests: leaves, duplications, losses, extra lineages,
  heights, lengths and cherries. All tests pass; for example, the mean gene-tree total length
  is 14.61 in fast mode and 14.62 in SimPhy.
  `tests/unif_test.c` checks the fast-mode probabilities against 2048-bit MPFR; the worst
  relative error is 2.6e-10.

To reproduce the checks:

```bash
git submodule update --init      # fetches SimPhy v1.0.2 into upstream/SimPhy
make oracle                      # builds the reference binary upstream/SimPhy/bin/simphy_ref
make check                       # numerics unit test + byte-exact matrix (1 and 4 threads)
```

## Documentation

* [DOCS/01-reference-analysis.md](DOCS/01-reference-analysis.md): how SimPhy 1.0.2 works,
  where its time and memory go, and the defects found
* [DOCS/02-design.md](DOCS/02-design.md): design of SimPhy-X and its two modes
* [DOCS/03-usage.md](DOCS/03-usage.md): detailed usage, input formats, and behaviour where
  SimPhy fails
* [DOCS/04-verification.md](DOCS/04-verification.md): tests and benchmarks
* [DOCS/cookbook.sh](DOCS/cookbook.sh): annotated example commands

Known limitation: verbosity levels 5–6 (SimPhy's per-coalescence debugging output) print the
level-4 text only.

## Relationship to SimPhy

SimPhy-X is derived from SimPhy 1.0.2: tag `v1.0.2` of
[adamallo/SimPhy](https://github.com/adamallo/SimPhy), which is included as the git submodule
`upstream/SimPhy`. The simulator was rewritten in C, following SimPhy's source and its order
of random-number draws, so that compat mode reproduces SimPhy exactly. Around that core,
SimPhy-X adds:

* array-based trees,
* a thread pool,
* the fast-mode numerics, and
* deterministic handling of inputs on which SimPhy 1.0.2 crashes.

[NOTICE](NOTICE) lists the changes. Bug reports about SimPhy-X belong in this repository, not
upstream.

## Citation

If you use SimPhy-X, please cite both SimPhy-X and the original SimPhy paper. GitHub's
"Cite this repository" button reads [CITATION.cff](CITATION.cff).

```bibtex
@article{saha2026simphyx,
  author  = {Saha, Anik and Rahman, Ashrafur and Bayzid, Md. Shamsuzzoha},
  title   = {{SimPhy-X}: Scalable, parallel and {SimPhy}-compatible simulation of species, locus and gene trees},
  note    = {In preparation. Software: https://github.com/aaniksahaa/simphy-x},
  year    = {2026}
}

@article{mallo2016simphy,
  author  = {Mallo, Diego and De Oliveira Martins, Leonardo and Posada, David},
  title   = {{SimPhy}: Phylogenomic Simulation of Gene, Locus, and Species Trees},
  journal = {Systematic Biology},
  volume  = {65},
  number  = {2},
  pages   = {334--344},
  year    = {2016},
  doi     = {10.1093/sysbio/syv082}
}
```

## License

SimPhy-X is free software. You can redistribute it and modify it under the terms of the
[GNU General Public License](LICENSE), version 2 or, at your option, any later version. This is
the same license family as SimPhy. Copyright © 2015 Diego Mallo, Leonardo de Oliveira Martins
and David Posada (SimPhy). Copyright © 2026 Anik Saha (SimPhy-X modifications). See [NOTICE](NOTICE) for
third-party libraries in the static binary.
