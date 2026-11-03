/* SPDX-License-Identifier: GPL-3.0-only */
/* Gene tree simulation along a locus tree: multispecies coalescent (MSC) and
 * bounded multilocus coalescent (MLC).  Active gene lineages live in a single
 * pool of n_gleaves ints used as a stack of per-locus-node segments.
 *
 * The MLC simulation can run in three ways with identical random draws:
 *   inline  : everything (including the MPFR fallback) is computed at once;
 *   record  : all draws are made, the pair choices and sampled waiting times
 *             are stored in an sx_mlc_trace, and MPFR chains are deferred
 *             (recorded in an sx_chain_list, their times left NaN);
 *   replay  : no draws; the gene tree is rebuilt from the trace once the
 *             chains have been evaluated (bit-identical arithmetic order). */
#ifndef SX_GTREE_H
#define SX_GTREE_H
#include "sx_tree.h"
#include "sx_rng.h"
#include "sx_coal.h"

typedef struct {
    int *pool;       /* n_gleaves ids */
    int *seg_stack;  /* stack of segment offsets pushed by processed nodes */
    int cap_nodes, cap_pool;
} sx_lineage_pool;

void sx_pool_init(sx_lineage_pool *p);
void sx_pool_free(sx_lineage_pool *p);
void sx_pool_reserve(sx_lineage_pool *p, int n_lnodes, int n_gleaves);

typedef struct {
    int n_lnodes;
    int *n_coals;        /* per locus node: number of coalescences performed */
    int *ngen_off;       /* per locus node: offset into ngens */
    double *ngens;       /* sampled waiting times, flat */
    int n_ngens, cap_ngens;
    int *picks;          /* two ints per coalescence, in processing order */
    int n_picks, cap_picks;
    int cap_l;
} sx_mlc_trace;

void sx_trace_init(sx_mlc_trace *t);
void sx_trace_free(sx_mlc_trace *t);

typedef struct {
    sx_mpfr_chain_job job;   /* u/out pointers are resolved at run time */
    int lnode;
    int u_off;               /* into sx_chain_list::ubuf */
    int out_off;             /* into sx_mlc_trace::ngens */
} sx_chain_rec;

typedef struct {
    sx_chain_rec *v;
    int n, cap;
    double *ubuf; int nu, ucap;
} sx_chain_list;

void sx_chains_init(sx_chain_list *c);
void sx_chains_free(sx_chain_list *c);
void sx_chains_clear(sx_chain_list *c);
/* Evaluate every recorded chain (exact MPFR arithmetic), writing the times into trace->ngens. */
void sx_chains_run_all(sx_chain_list *c, sx_mlc_trace *trace);
/* Evaluate one chain (for finer-grained parallelism). */
void sx_chain_run_one(sx_chain_list *c, sx_mlc_trace *trace, int idx);

/* SimMSCGTree port. simlosses: number of g leaves per lost l leaf (0/1). */
long sx_sim_msc_gtree(sx_ltree *lt, sx_gtree *gt, sx_lineage_pool *pool, gsl_rng *seed, int *n_lcoals, int simlosses, int verbosity, int collapse);

/* SimMLCGTree port (compat numerics).
 * trace == NULL: inline mode; the gene tree is complete on return (collapsed if collapse).
 * trace != NULL: record mode; chains may be NULL (MPFR inline) or not (deferred).
 *   The tree is NOT collapsed; call sx_sim_mlc_replay (if chains were
 *   recorded) and then sx_gtree_collapse yourself. */
long sx_sim_mlc_gtree(sx_ltree *lt, sx_gtree *gt, sx_lineage_pool *pool, float epsilon_brent, gsl_rng *seed, int *tn_lcoals, int verbosity, int collapse,
                      sx_chain_list *chains, sx_mlc_trace *trace);
long sx_sim_mlc_replay(sx_ltree *lt, sx_gtree *gt, sx_lineage_pool *pool, const sx_mlc_trace *trace);
/* Fast mode: same model, bounded waiting times sampled by uniformisation (sx_unif). */
struct sx_unif_ws;
long sx_sim_mlc_gtree_fast(sx_ltree *lt, sx_gtree *gt, sx_lineage_pool *pool, float epsilon_brent, gsl_rng *seed, int *tn_lcoals, int collapse, struct sx_unif_ws *ws, int compat_probs);

/* Rateheter_GTbranchspec: one gamma per node in array order. */
void sx_gtree_rateheter(sx_gtree *gt, double alpha, gsl_rng *seed);
/* Draw the multipliers only (array order), to be applied later with sx_gtree_apply_multipliers. */
void sx_gtree_draw_multipliers(int n_nodes, double alpha, gsl_rng *seed, double *out);
void sx_gtree_apply_multipliers(sx_gtree *gt, const double *m);

#endif
