# simphy-x — usage

## Build

```
make            # → bin/simphy-x  (GSL, MPFR, GMP, SQLite linked statically from deps/)
make clean
```

Requirements: gcc, GSL development headers (`libgsl-dev`), and the MPFR/GMP/SQLite headers and
static libraries. On this machine the latter were extracted from the Ubuntu packages into
`deps/` (`apt-get download libmpfr-dev libgmp-dev libsqlite3-dev` then `dpkg -x`), no root
needed; point `DEPS=` at another prefix if required.

## Command line

`simphy-x` accepts SimPhy 1.0.2's command line unchanged (see `simphy-x -h` for the full list:
`-RS -RL -RG`, `-S -SB -SD -ST -SL -SO -SI -SP -SU -SG`, `-LB -LD -LT -LG -LK -LL -LS`,
`-GB -GD -GT -GG -GP`, `-HS -HL -HH -HG`, `-CS -CE`, `-V`, `-O -OT -OM -OD -OP -OC -OL -ON -OS`)
and produces the same directory layout (`<prefix>/<replicate>/s_tree.trees`, `l_trees.trees`,
`g_treesNNN.trees`, `*.mapsl`, `*.maplg`, `bounded_locus_subtrees.out`, `stats.txt`,
`<prefix>.params`, `<prefix>.command`, `<prefix>.db`).

New options:

| option | meaning |
|---|---|
| `-XM compat` (default) | byte-identical output to SimPhy 1.0.2 for the same command and seed, at any thread count |
| `-XM fast` | same model and distributions; independent random stream per (replicate, locus); no arbitrary-precision arithmetic in the common path; fully parallel |
| `-XT n` | number of threads (default: all cores). Output never depends on `n` |
| `-XV`, `--version` | print the version |

Long forms `--mode`, `--threads`, `--help` are also accepted.

Environment variables keep SimPhy's meaning: `SIMPHY_MAXLEAVES` (cap on locus-tree lineages
per species branch and on expected locus-tree leaves; raise it for large duplication rates),
`SIMPHY_MAXIT`, `SIMPHY_MAXNAME`, `SIMPHY_FLOATPREC`. Diagnostics: `SIMPHY_X_DEBUG=1` prints
fast-mode sampler counters; `SIMPHY_X_PROBLOG=file` logs compat-mode lineage-count
probability evaluations.

Not supported yet: Nexus tree files (`-SR`, `-LR`), fixed locus trees (`-L`), configuration
files (`-I`). The program exits with a settings error for these.

## Which mode to use

* Reproducing an existing dataset, or when byte-for-byte agreement with SimPhy matters:
  `compat`. The multispecies coalescent pipelines gain mostly memory (O(n) instead of O(n²));
  the duplication/loss pipelines gain the parallel evaluation of the arbitrary-precision
  chains (4–5× on 16 cores; the remaining part is inherently sequential because SimPhy uses
  one random stream).
* New datasets, large taxon counts, duplication/loss models: `fast`. Same distributions
  (verified statistically, see 04-verification.md), 50–150× faster than SimPhy on
  duplication/loss scenarios, deterministic for a given seed regardless of thread count.
  Fast-mode output is *not* the same sample path as SimPhy's for the same seed.

## Pipeline integration

`ref/normal-phylo-sim/run_simulator.sh` and `ref/duploss-simulator/duploss-sim.sh` can use
`simphy-x` as a drop-in replacement for `simphy_lnx64` (same flags, same output files). To use
fast mode append `-XM fast` to the SimPhy invocation; `-XT` defaults to all cores.

## Seeds and determinism

* compat: the seed is passed to GSL's mt19937 exactly as SimPhy does; identical output.
* fast: the seed and the (replicate, locus) indices are hashed into independent 64-bit
  stream seeds (splitmix64) for xoshiro256** generators; the species tree of replicate r uses
  stream (seed, r), locus l of replicate r uses stream (seed, r, l). Changing the number of
  loci or threads does not change the species trees or the other loci.
