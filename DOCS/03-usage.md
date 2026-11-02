# simphy-x — usage

## Build

```
make            # → bin/simphy-x          fully static, portable (the binary to distribute)
make dynamic    # → bin/simphy-x-dynamic  linked against the host glibc/libm
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
| `-XR 1` | compat only: where SimPhy 1.0.2 aborts with "MPFR out of precision" (bounded subtrees with several hundred lineages), retry with doubled precision instead. Output is identical up to that locus and deterministic afterwards, but no longer comparable to SimPhy, which has no output there. A note is printed on stderr when it triggers |
| `-XV`, `--version` | print the version |

Long forms `--mode`, `--threads`, `--help` are also accepted.

Environment variables keep SimPhy's meaning: `SIMPHY_MAXLEAVES` (cap on locus-tree lineages
per species branch and on expected locus-tree leaves; raise it for large duplication rates),
`SIMPHY_MAXIT`, `SIMPHY_MAXNAME`, `SIMPHY_FLOATPREC`. Diagnostics: `SIMPHY_X_DEBUG=1` prints
fast-mode sampler counters; `SIMPHY_X_PROBLOG=file` logs compat-mode lineage-count
probability evaluations.

Configuration files (`-I`) are supported (one `-CODE value` per line, `#`/`//` comments;
command-line options take precedence; the file is copied to `<prefix>.conf` when `-OC 1`).

## Fixed trees as input

All four SimPhy input forms are implemented and verified byte-exact against the reference.

| Option | Input | Notes |
|---|---|---|
| `-S newick` | one fixed species tree, extended Newick | per-node codes `#Ne`, `/n_individuals`, `*mu_mult`, `~gtime_mult`; the tree is used for every replicate |
| `-SR file.nex` | Nexus file with one species tree per replicate | `[&pop_size=,n_ind=,u_mult=,g_mult=]` comments after a node; `-RS` is set to the number of trees in the file |
| `-L newick` | one fixed locus tree, extended Newick | leaf `Name_paralog`, node kind `%k` (1 duplication, 2 loss, 3 transfer, 4 transfer receptor, 5 conversion, 6 conversion receptor), `/n_individuals`, `#Ne`, `*mu_mult`, `~gtime_mult`; no species tree exists, so no `s_tree.trees` and no mappings |
| `-LR file.nex` | Nexus file with one locus tree per locus | `[&pop_size=,kind_n=,paralog=,n_ind=,u_mult=,g_mult=]`; `-RL` is set to the number of trees in the file |

Nexus files must start with a line `#NEXUS`, contain a line `BEGIN TREES;` (any case), one
`tree name = (...);` line per tree and an `End;` line (see `tests/inputs/*.nex`). Fixed
locus trees are parsed once per locus and the lineage-count probabilities are computed with
the compat arithmetic, so a fixed locus tree with duplications gives the same gene trees as
SimPhy.

Where SimPhy 1.0.2 itself fails on these inputs simphy-x keeps going (there is nothing to be
byte-identical to):

* `-L`/`-LR` with `-HL` (gene-family rate heterogeneity) crashes SimPhy (it multiplies the
  species tree's rate, and there is no species tree); simphy-x uses `-SU` as the base rate.
* `-L`/`-LR` with `-OT 1` (branches in time units) crashes SimPhy; simphy-x writes the locus
  and gene trees in time units from the locus tree's generation-time multipliers.
* `-L`/`-LR` with `-OM 1`: SimPhy prints the "mapping outputs will not be generated" warning
  but, because of a pointer/value mix-up, leaves mapping on; simphy-x prints the same warning
  and turns mapping off.
* `-SR`/`-LR`: SimPhy's Nexus line reader appends to an uninitialised buffer, so the tree
  string it parses depends on heap contents (it comes out corrupted for some verbosity levels
  and file layouts). simphy-x reads the line deterministically. It matches SimPhy whenever
  SimPhy's read happens to be clean, which is the case on the verification matrix.
* `-S` together with `-SR` (or `-L` with `-LR`): SimPhy frees the Newick string and then keeps
  using it; simphy-x prints SimPhy's warning and uses the file.

## Program text (stdout)

After the simphy-x banner (ASCII art, version, mode, thread count, separator line) the
program text is SimPhy's, character for character, for verbosity 0–4: the settings block,
replicate/locus/gene-tree parameter blocks, progress messages, and the trees echoed at
verbosity 4. Loci finish out of order in parallel runs, so their text is buffered and printed
in locus order. Verbosity 5–6 (SimPhy's per-coalescence debugging output) prints the level-4
text only. `tests/verify_compat.sh` checks stdout too when `SX_CHECK_STDOUT=1`.

## Portable binary

`make` builds `bin/simphy-x` as a fully static x86-64 Linux executable (glibc, GSL,
MPFR, GMP, SQLite all linked in, 3.5 MB). It runs on any x86-64 Linux without dependencies and
gives byte-identical output on every machine, because the `exp`/`log` implementation is baked
in (`simphy-x-dynamic` and SimPhy itself use the host's libm, which can differ between glibc
versions). It needs no shared libraries at all: `ldd` reports "not a dynamic executable". The
suppressed link warning concerns SQLite's `dlopen` for loadable extensions, which are not used.
It is x86-64 Linux only; other platforms build from source with `make`.

## Which mode to use

* Reproducing an existing dataset, or when byte-for-byte agreement with SimPhy matters:
  `compat`. The multispecies coalescent pipelines gain mostly memory (O(n) instead of O(n²));
  the duplication/loss pipelines gain the parallel evaluation of the arbitrary-precision
  chains (4–5× on 16 cores; the remaining part is inherently sequential because SimPhy uses
  one random stream).
* At large scale SimPhy 1.0.2 cannot run some duplication/loss cases at all: with 2000 taxa
  its 512-bit MPFR sampler aborts ("MPFR out of precision"), and compat mode aborts identically
  after the same 495 MB → 33 MB of memory; fast mode completes them (2001-leaf gene trees,
  7 MB).
* New datasets, large taxon counts, duplication/loss models: `fast`. Same distributions
  (verified statistically, see 04-verification.md), 50–150× faster than SimPhy on
  duplication/loss scenarios, deterministic for a given seed regardless of thread count.
  Fast-mode output is *not* the same sample path as SimPhy's for the same seed.

## Pipeline integration

Scripts that call SimPhy 1.0.2 (`simphy_lnx64`) can use `simphy-x` as a drop-in replacement
(same flags, same output files). To use
fast mode append `-XM fast` to the SimPhy invocation; `-XT` defaults to all cores.

## Seeds and determinism

* compat: the seed is passed to GSL's mt19937 exactly as SimPhy does; identical output.
* fast: the seed and the (replicate, locus) indices are hashed into independent 64-bit
  stream seeds (splitmix64) for xoshiro256** generators; the species tree of replicate r uses
  stream (seed, r), locus l of replicate r uses stream (seed, r, l). Changing the number of
  loci or threads does not change the species trees or the other loci.
