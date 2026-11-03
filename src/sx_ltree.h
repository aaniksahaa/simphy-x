/* SPDX-License-Identifier: GPL-3.0-only */
/* Locus tree simulation (SimBDLHTree port) and the lineage-count probability
 * machinery used by the bounded multilocus coalescent. */
#ifndef SX_LTREE_H
#define SX_LTREE_H
#include "sx_tree.h"
#include "sx_rng.h"

typedef struct {
    int losses, dups, transfr, gc, leaves, gleaves, ntrials;
} sx_ltree_stats;

/* node_ptrs: scratch of MAX_LEAVES ints owned by the caller. */
long sx_sim_bdlh_tree(sx_stree *st, sx_ltree **lt, int *node_ptrs, double b_rate, double d_rate, double h_rate, double gc_rate, int t_kind,
                      gsl_rng *seed, int min_lleaves, int min_lsleaves, int verbosity, sx_ltree_stats *stats, int calc_probs);

struct sx_tpool;
/* Optional thread pool for the compat probability pass (rows are computed in parallel with
 * unchanged per-row arithmetic; results are bit-identical). NULL = sequential. */
void sx_calc_probs_set_pool(struct sx_tpool *pool);
/* CalcProbsNLineagesLTree (compat arithmetic). Returns fmax_nlin of node. */
int sx_calc_probs_nlineages(sx_ltree *t, int node, int Ne, int in_subtree, int verbosity);
/* SampleNLineagesCDFLTree (compat). */
void sx_sample_nlineages_cdf(sx_ltree *t, int node, int n_gleaves, int Ne, int verbosity, gsl_rng *seed);

double sx_expected_pruned_ltree_nleaves(const sx_stree *st, double b_rate, double d_rate);

#endif

/* ---- fast mode (uniformisation) ---- */
struct sx_unif_ws;
/* Same semantics as sx_calc_probs_nlineages but all probabilities computed by
 * uniformisation (no cancellation).  Returns fmax_nlin. */
int sx_calc_probs_nlineages_fast(sx_ltree *t, int node, int Ne, int in_subtree, struct sx_unif_ws *ws);
/* Same semantics as sx_sample_nlineages_cdf (draw order: one uniform for kt,
 * one for the split) with uniformisation probabilities. */
void sx_sample_nlineages_fast(sx_ltree *t, int node, int Ne, gsl_rng *seed, struct sx_unif_ws *ws);
