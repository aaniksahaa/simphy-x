/* SPDX-License-Identifier: GPL-3.0-only */
#include "sx_gtree.h"
#include "sx_ltree.h"
#include "sx_unif.h"

void sx_pool_init(sx_lineage_pool *p) { memset(p, 0, sizeof(*p)); }
void sx_pool_free(sx_lineage_pool *p) { free(p->pool); free(p->seg_stack); memset(p, 0, sizeof(*p)); }
void sx_pool_reserve(sx_lineage_pool *p, int n_lnodes, int n_gleaves)
{
    if (n_lnodes + 1 > p->cap_nodes) {
        p->seg_stack = sx_realloc(p->seg_stack, sizeof(int) * (size_t)(n_lnodes + 1));
        p->cap_nodes = n_lnodes + 1;
    }
    if (n_gleaves + 1 > p->cap_pool) {
        p->pool = sx_realloc(p->pool, sizeof(int) * (size_t)(n_gleaves + 1));
        p->cap_pool = n_gleaves + 1;
    }
}

void sx_trace_init(sx_mlc_trace *t) { memset(t, 0, sizeof(*t)); }
void sx_trace_free(sx_mlc_trace *t) { free(t->n_coals); free(t->ngen_off); free(t->ngens); free(t->picks); memset(t, 0, sizeof(*t)); }
static void trace_begin(sx_mlc_trace *t, int n_lnodes)
{
    if (n_lnodes > t->cap_l) {
        t->n_coals = sx_realloc(t->n_coals, sizeof(int) * (size_t)n_lnodes);
        t->ngen_off = sx_realloc(t->ngen_off, sizeof(int) * (size_t)n_lnodes);
        t->cap_l = n_lnodes;
    }
    t->n_lnodes = n_lnodes; t->n_ngens = 0; t->n_picks = 0;
}
static int trace_push_ngens(sx_mlc_trace *t, const double *v, int n)
{
    if (t->n_ngens + n > t->cap_ngens) {
        int nc = t->cap_ngens ? t->cap_ngens : 256;
        while (nc < t->n_ngens + n) nc *= 2;
        t->ngens = sx_realloc(t->ngens, sizeof(double) * (size_t)nc);
        t->cap_ngens = nc;
    }
    int off = t->n_ngens;
    if (n > 0) memcpy(t->ngens + off, v, sizeof(double) * (size_t)n);
    t->n_ngens += n;
    return off;
}
static inline void trace_push_pick(sx_mlc_trace *t, int a, int b)
{
    if (t->n_picks + 2 > t->cap_picks) {
        int nc = t->cap_picks ? t->cap_picks : 256;
        while (nc < t->n_picks + 2) nc *= 2;
        t->picks = sx_realloc(t->picks, sizeof(int) * (size_t)nc);
        t->cap_picks = nc;
    }
    t->picks[t->n_picks++] = a; t->picks[t->n_picks++] = b;
}

void sx_chains_init(sx_chain_list *c) { memset(c, 0, sizeof(*c)); }
void sx_chains_free(sx_chain_list *c) { free(c->v); free(c->ubuf); memset(c, 0, sizeof(*c)); }
void sx_chains_clear(sx_chain_list *c) { c->n = 0; c->nu = 0; }

static sx_chain_rec *chains_add(sx_chain_list *c, int n_coals)
{
    if (c->n == c->cap) { c->cap = c->cap ? c->cap * 2 : 64; c->v = sx_realloc(c->v, sizeof(sx_chain_rec) * (size_t)c->cap); }
    if (c->nu + n_coals > c->ucap) {
        int nc = c->ucap ? c->ucap : 1024;
        while (nc < c->nu + n_coals) nc *= 2;
        c->ubuf = sx_realloc(c->ubuf, sizeof(double) * (size_t)nc);
        c->ucap = nc;
    }
    sx_chain_rec *r = &c->v[c->n++];
    memset(r, 0, sizeof(*r));
    r->job.n_coals = n_coals;
    r->u_off = c->nu;
    c->nu += n_coals;
    return r;
}

void sx_chain_run_one(sx_chain_list *c, sx_mlc_trace *trace, int idx)
{
    sx_chain_rec *r = &c->v[idx];
    r->job.u = c->ubuf + r->u_off;
    r->job.out = trace->ngens + r->out_off;
    sx_mpfr_chain_run(&r->job);
}
void sx_chains_run_all(sx_chain_list *c, sx_mlc_trace *trace)
{
    for (int i = 0; i < c->n; ++i) sx_chain_run_one(c, trace, i);
}

/* ---- pick source: RNG (optionally recording) or replay ---- */
typedef struct { gsl_rng *seed; sx_mlc_trace *rec; const sx_mlc_trace *rep; int cursor; } pick_src;

static inline void get_picks(pick_src *ps, int avail_before, int *p1, int *p2)
{
    if (ps->rep) {
        *p1 = ps->rep->picks[ps->cursor++];
        *p2 = ps->rep->picks[ps->cursor++];
        return;
    }
    *p1 = (int)(gsl_rng_uniform(ps->seed) * avail_before);
    int avail_after = avail_before - 1;
    if (avail_after > 1) *p2 = (int)(gsl_rng_uniform(ps->seed) * avail_after);
    else *p2 = 0;
    if (ps->rec) trace_push_pick(ps->rec, *p1, *p2);
}

/* Common coalescence step: mirrors SimPhy exactly. */
static inline void coalesce_pair(sx_gtree *gt, sx_ltree *lt, int w_lnode, int *seg, int *avail_leaves, int *next_avail_inode,
                                 double current_ngen, double p_mu, pick_src *ps, int *anc_out)
{
    int p1, p2;
    get_picks(ps, *avail_leaves, &p1, &p2);
    int off1 = seg[p1];
    --(*avail_leaves);
    seg[p1] = seg[*avail_leaves];
    int off2 = seg[p2];
    int anc = *next_avail_inode;
    seg[p2] = anc;
    if (gt->sp_index[off1] == gt->sp_index[off2]) gt->sp_index[anc] = gt->sp_index[off1];
    ++(*next_avail_inode);
    --lt->n_lin[w_lnode];
    SX_CHILD(gt, anc, gt->n_child[anc]) = off1;
    SX_CHILD(gt, anc, gt->n_child[anc] + 1) = off2;
    gt->anc[off1] = anc; gt->anc[off2] = anc;
    gt->contl[anc] = w_lnode; gt->conts[anc] = lt->conts[w_lnode];
    gt->n_child[anc] += 2;
    gt->n_gen[anc] = current_ngen;
    double top = lt->n_gen[w_lnode];
    gt->gen_length[off1] = gt->n_gen[off1] - gt->n_gen[anc];
    gt->bl[off1] += ((gt->n_gen[off1] > top ? top : gt->n_gen[off1]) - gt->n_gen[anc]) * p_mu;
    gt->gen_length[off2] = gt->n_gen[off2] - gt->n_gen[anc];
    gt->bl[off2] += ((gt->n_gen[off2] > top ? top : gt->n_gen[off2]) - gt->n_gen[anc]) * p_mu;
    gt->paralog[anc] = lt->paralog[w_lnode];
    *anc_out = anc;
}

static inline void hand_to_parent(sx_gtree *gt, sx_ltree *lt, int w_lnode, const int *seg, int n_left, double min_ngen, double p_mu)
{
    for (int j = 0; j < n_left; ++j) {
        int g = seg[j];
        gt->bl[g] += (gt->n_gen[g] > lt->n_gen[w_lnode] ? lt->gen_length[w_lnode] : gt->n_gen[g] - min_ngen) * p_mu;
    }
}

static inline void place_leaves(sx_gtree *gt, sx_ltree *lt, int w_lnode, int *seg, int n, int *next_leaf)
{
    for (int j = 0; j < n; ++j) {
        int g = (*next_leaf)++;
        seg[j] = g;
        gt->sp_index[g] = lt->sp_index[w_lnode];
        gt->paralog[g] = lt->paralog[w_lnode];
        gt->replica[g] = j;
        gt->contl[g] = w_lnode;
        gt->conts[g] = lt->conts[w_lnode];
        gt->n_gen[g] = lt->n_gen[w_lnode];
    }
}

/* ======================= MSC ======================= */
long sx_sim_msc_gtree(sx_ltree *lt, sx_gtree *gt, sx_lineage_pool *pool, gsl_rng *seed, int *n_lcoals, int simlosses, int verbosity, int collapse)
{
    (void)verbosity;
    if (simlosses > 1) return SX_UNEXPECTED_VALUE;
    sx_gtree_reset(gt);
    sx_pool_reserve(pool, lt->n_nodes, lt->n_gleaves);
    for (int i = 0; i < lt->n_nodes; ++i) {
        if (lt->n_child[i] == 0) {
            switch (lt->kind[i]) {
            case LOSS: case RTRFR: case RGC: lt->n_lin[i] = simlosses; break;
            default: lt->n_lin[i] = lt->conts[i] >= 0 ? lt->species_tree->n_replicas[lt->conts[i]] : lt->n_ilin[i]; break;
            }
        } else lt->n_lin[i] = 0;
    }
    gt->species_tree = lt->species_tree; gt->locus_tree = lt; gt->gen_time = lt->gen_time;
    pick_src ps = { seed, NULL, NULL, 0 };

    int next_leaf = 0, next_avail_inode = lt->n_gleaves, sp = 0, top = 0, anc_gnode = -1;
    *n_lcoals = 0;
    for (int i = 0; i < lt->n_nodes; ++i) {
        int w_lnode = i, nch = lt->n_child[i];
        int seg_off;
        if (nch == 0) { seg_off = top; place_leaves(gt, lt, w_lnode, pool->pool + seg_off, lt->n_lin[w_lnode], &next_leaf); }
        else { sp -= nch; seg_off = pool->seg_stack[sp]; }
        int *seg = pool->pool + seg_off;
        int k = lt->n_lin[w_lnode];
        if (k == 0) { pool->seg_stack[sp++] = seg_off; top = seg_off; continue; }
        int anc_lnode = lt->anc[w_lnode];
        int avail_leaves = k;
        int p_Ne = lt->Ne_b[w_lnode] != 0 ? lt->Ne_b[w_lnode] : lt->Ne;
        double p_mu = lt->mu * lt->mu_mult[w_lnode];
        double current_ngen = lt->n_gen[w_lnode];
        double min_ngen = current_ngen - lt->gen_length[w_lnode];
        int max_coals = k - 1;
        for (int j = 0; j < max_coals; ++j) {
            double sampled_ngen = ((-log(gsl_rng_uniform_pos(seed))) * (2 * p_Ne)) / (k * (k - 1));
            if (sampled_ngen <= 0) return SX_UNEXPECTED_VALUE;
            current_ngen -= sampled_ngen;
            --k;
            if (current_ngen < min_ngen && anc_lnode != -1) { (*n_lcoals) += (avail_leaves - 1); break; }
            coalesce_pair(gt, lt, w_lnode, seg, &avail_leaves, &next_avail_inode, current_ngen, p_mu, &ps, &anc_gnode);
        }
        int n_left = lt->n_lin[w_lnode];
        if (anc_lnode != -1) { hand_to_parent(gt, lt, w_lnode, seg, n_left, min_ngen, p_mu); lt->n_lin[anc_lnode] += n_left; }
        else gt->root = anc_gnode;
        pool->seg_stack[sp++] = seg_off;
        top = seg_off + n_left;
    }
    if (collapse == 1) sx_gtree_collapse(gt, 1);
    return SX_NO_ERROR;
}

/* ======================= MLC ======================= */
/* Shared per-node loop. mode: 0 inline, 1 record, 2 replay. */
static long mlc_loop(sx_ltree *lt, sx_gtree *gt, sx_lineage_pool *pool, float epsilon_brent, gsl_rng *seed, int *tn_lcoals, int verbosity,
                     sx_chain_list *chains, sx_mlc_trace *rec, const sx_mlc_trace *rep, sx_unif_ws *fast)
{
    int mode = rep ? 2 : (rec ? 1 : 0);
    sx_gtree_reset(gt);
    sx_pool_reserve(pool, lt->n_nodes, lt->n_gleaves);
    /* Active lineages: leaves start with their gene leaves, internal nodes
     * with none.  n_ilin of bounded internal nodes (set by the count sampling)
     * must be preserved: the MPFR fallback reads it. */
    for (int i = 0; i < lt->n_nodes; ++i) lt->n_lin[i] = lt->n_child[i] == 0 ? lt->n_ilin[i] : 0;
    gt->species_tree = lt->species_tree; gt->locus_tree = lt; gt->gen_time = lt->gen_time;
    if (rec) trace_begin(rec, lt->n_nodes);
    pick_src ps = { seed, rec, rep, 0 };

    double *sampled_ngens = sx_calloc((size_t)lt->n_gleaves + 1, sizeof(double));
    int next_leaf = 0, next_avail_inode = lt->n_gleaves, sp = 0, top = 0, anc_gnode = -1;
    int tn = 0;

    for (int i = 0; i < lt->n_nodes; ++i) {
        int w_lnode = i, nch = lt->n_child[i];
        int seg_off;
        if (nch == 0) { seg_off = top; place_leaves(gt, lt, w_lnode, pool->pool + seg_off, lt->n_lin[w_lnode], &next_leaf); }
        else { sp -= nch; seg_off = pool->seg_stack[sp]; }
        int *seg = pool->pool + seg_off;
        int n_lin = lt->n_lin[w_lnode];
        if (n_lin == 0) {
            if (rec) { rec->n_coals[i] = 0; rec->ngen_off[i] = rec->n_ngens; }
            pool->seg_stack[sp++] = seg_off; top = seg_off; continue;
        }
        int n_coals, no_counts;
        if (lt->n_olin[w_lnode] == 0) { n_coals = n_lin - 1; lt->n_ilin[w_lnode] = n_lin; no_counts = 1; }
        else { no_counts = 0; n_coals = n_lin - lt->n_olin[w_lnode]; tn += lt->n_olin[w_lnode] - 1; }

        int anc_lnode = lt->anc[w_lnode];
        int avail_leaves = n_lin;
        int p_Ne = lt->Ne_b[w_lnode] != 0 ? lt->Ne_b[w_lnode] : lt->Ne;
        double p_mu = lt->mu * lt->mu_mult[w_lnode];
        double current_ngen = lt->n_gen[w_lnode];
        double min_ngen = current_ngen - lt->gen_length[w_lnode];

        if (mode == 2) {
            n_coals = rep->n_coals[i];
            if (n_coals > 0) memcpy(sampled_ngens, rep->ngens + rep->ngen_off[i], sizeof(double) * (size_t)n_coals);
        } else if (n_coals > 0) {
            int deferred_here = 0;
            if (lt->gen_length[w_lnode] != 0) {
                int j;
                if (no_counts == 1) {
                    j = 0;
                    sampled_ngens[0] = ((-log(gsl_rng_uniform_pos(seed))) * (2 * p_Ne)) / (avail_leaves * (avail_leaves - 1));
                    while (j < n_coals && current_ngen - sampled_ngens[j] > min_ngen) {
                        current_ngen -= sampled_ngens[j];
                        --avail_leaves;
                        ++j;
                        sampled_ngens[j] = ((-log(gsl_rng_uniform_pos(seed))) * (2 * p_Ne)) / (avail_leaves * (avail_leaves - 1));
                    }
                    n_coals = j;
                } else if (fast) {
                    long frc = sx_sample_bounded_waits(fast, avail_leaves, lt->n_olin[w_lnode], current_ngen - min_ngen, (double)p_Ne, epsilon_brent, seed, sampled_ngens);
                    if (frc != SX_NO_ERROR) { free(sampled_ngens); return frc; }
                } else {
                    for (j = 0; j < n_coals; ++j) {
                        sampled_ngens[j] = sx_sample_coal_time_mlc(avail_leaves, lt->n_olin[w_lnode], current_ngen - min_ngen, p_Ne, epsilon_brent, gsl_rng_uniform_pos(seed), verbosity);
                        current_ngen -= sampled_ngens[j];
                        if (current_ngen - epsilon_brent <= min_ngen) break;
                        --avail_leaves;
                    }
                    if (sampled_ngens[n_coals - 1] <= 0 || current_ngen - epsilon_brent <= min_ngen) {
                        current_ngen = lt->n_gen[w_lnode];
                        avail_leaves = lt->n_ilin[w_lnode];
                        if (chains && rec) {
                            sx_chain_rec *cr = chains_add(chains, n_coals);
                            cr->lnode = w_lnode;
                            cr->job.n_ilin = avail_leaves; cr->job.n_olin = lt->n_olin[w_lnode]; cr->job.pop_size = p_Ne;
                            cr->job.n_gen = lt->n_gen[w_lnode]; cr->job.min_ngen = min_ngen; cr->job.epsilon = epsilon_brent; cr->job.verbosity = verbosity;
                            double *u = chains->ubuf + cr->u_off;
                            for (j = 0; j < n_coals; ++j) u[j] = gsl_rng_uniform_pos(seed);
                            for (j = 0; j < n_coals; ++j) sampled_ngens[j] = NAN;
                            cr->out_off = rec->n_ngens;   /* this node's ngens are pushed next */
                            deferred_here = 1;
                        } else {
                            for (j = 0; j < n_coals; ++j) {
                                sampled_ngens[j] = sx_mpfr_sample_coal_time_mlc(avail_leaves, lt->n_olin[w_lnode], current_ngen - min_ngen, p_Ne, epsilon_brent, gsl_rng_uniform_pos(seed), verbosity, 512);
                                current_ngen -= sampled_ngens[j];
                                --avail_leaves;
                            }
                        }
                    }
                }
                if (!deferred_here && n_coals > 1 && (sampled_ngens[n_coals - 1] <= 0 || current_ngen - epsilon_brent <= min_ngen)) {
                    free(sampled_ngens);
                    return SX_UNEXPECTED_VALUE;
                }
            } else {
                for (int j = 0; j < n_coals; ++j) {
                    sampled_ngens[j] = ((-log(gsl_rng_uniform_pos(seed))) * (2 * p_Ne)) / (avail_leaves * (avail_leaves - 1));
                    current_ngen -= sampled_ngens[j];
                    --avail_leaves;
                }
            }
        }
        if (rec) { rec->n_coals[i] = n_coals; rec->ngen_off[i] = trace_push_ngens(rec, sampled_ngens, n_coals); }

        current_ngen = lt->n_gen[w_lnode];
        avail_leaves = n_lin;
        min_ngen = current_ngen - lt->gen_length[w_lnode];
        for (int j = 0; j < n_coals; ++j) {
            current_ngen -= sampled_ngens[j];
            coalesce_pair(gt, lt, w_lnode, seg, &avail_leaves, &next_avail_inode, current_ngen, p_mu, &ps, &anc_gnode);
        }
        if (no_counts == 1) tn += (avail_leaves - 1);

        int n_left = lt->n_lin[w_lnode];
        if (anc_lnode != -1) { hand_to_parent(gt, lt, w_lnode, seg, n_left, min_ngen, p_mu); lt->n_lin[anc_lnode] += n_left; }
        else gt->root = anc_gnode;
        pool->seg_stack[sp++] = seg_off;
        top = seg_off + n_left;
    }
    free(sampled_ngens);
    if (tn_lcoals) *tn_lcoals = tn;
    return SX_NO_ERROR;
}

long sx_sim_mlc_gtree(sx_ltree *lt, sx_gtree *gt, sx_lineage_pool *pool, float epsilon_brent, gsl_rng *seed, int *tn_lcoals, int verbosity, int collapse,
                      sx_chain_list *chains, sx_mlc_trace *trace)
{
    /* MatchTreesMLC: bounded roots get n_olin=1, everything else 0 */
    for (int i = 0; i < lt->n_nodes; ++i) {
        int anc = lt->anc[i];
        if (anc != -1 && ((lt->kind[anc] == DUP || lt->kind[anc] == TRFR || lt->kind[anc] == GC) && SX_CHILD(lt, anc, 1) == i) && lt->fmax_nlin[i] > 0)
            lt->n_olin[i] = 1;
        else
            lt->n_olin[i] = 0;
        if (lt->n_child[i] != 0) lt->n_ilin[i] = 0;
    }
    if (getenv("SIMPHY_X_COMPAT_EXACTPROBS")) {
        static __thread sx_unif_ws *dws = NULL;
        if (!dws) { dws = sx_malloc(sizeof(*dws)); sx_unif_ws_init(dws); }
        sx_sample_nlineages_fast(lt, lt->root, lt->Ne, seed, dws);
    } else
        sx_sample_nlineages_cdf(lt, lt->root, lt->n_gleaves, lt->Ne, verbosity, seed);
    long rc = mlc_loop(lt, gt, pool, epsilon_brent, seed, tn_lcoals, verbosity, trace ? chains : NULL, trace, NULL, NULL);
    if (rc != SX_NO_ERROR) return rc;
    if (!trace && collapse == 1) sx_gtree_collapse(gt, 1);
    return SX_NO_ERROR;
}

long sx_sim_mlc_replay(sx_ltree *lt, sx_gtree *gt, sx_lineage_pool *pool, const sx_mlc_trace *trace)
{
    return mlc_loop(lt, gt, pool, 0.0f, NULL, NULL, 0, NULL, NULL, trace, NULL);
}

long sx_sim_mlc_gtree_fast(sx_ltree *lt, sx_gtree *gt, sx_lineage_pool *pool, float epsilon_brent, gsl_rng *seed, int *tn_lcoals, int collapse, sx_unif_ws *ws, int compat_probs)
{
    for (int i = 0; i < lt->n_nodes; ++i) {
        int anc = lt->anc[i];
        if (anc != -1 && ((lt->kind[anc] == DUP || lt->kind[anc] == TRFR || lt->kind[anc] == GC) && SX_CHILD(lt, anc, 1) == i) && lt->fmax_nlin[i] > 0)
            lt->n_olin[i] = 1;
        else
            lt->n_olin[i] = 0;
        if (lt->n_child[i] != 0) lt->n_ilin[i] = 0;
    }
    if (compat_probs) sx_sample_nlineages_cdf(lt, lt->root, lt->n_gleaves, lt->Ne, 0, seed);
    else sx_sample_nlineages_fast(lt, lt->root, lt->Ne, seed, ws);
    long rc = mlc_loop(lt, gt, pool, epsilon_brent, seed, tn_lcoals, 0, NULL, NULL, NULL, ws);
    if (rc != SX_NO_ERROR) return rc;
    if (collapse == 1) sx_gtree_collapse(gt, 1);
    return SX_NO_ERROR;
}

void sx_gtree_rateheter(sx_gtree *gt, double alpha, gsl_rng *seed)
{
    for (int i = 0; i < gt->n_nodes; ++i) gt->bl[i] *= gsl_ran_gamma(seed, alpha, 1 / alpha);
}
void sx_gtree_draw_multipliers(int n_nodes, double alpha, gsl_rng *seed, double *out)
{
    for (int i = 0; i < n_nodes; ++i) out[i] = gsl_ran_gamma(seed, alpha, 1 / alpha);
}
void sx_gtree_apply_multipliers(sx_gtree *gt, const double *m)
{
    for (int i = 0; i < gt->n_nodes; ++i) gt->bl[i] *= m[i];
}
