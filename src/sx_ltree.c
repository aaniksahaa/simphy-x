/* SPDX-License-Identifier: GPL-3.0-only */
#include "sx_ltree.h"
#include "sx_coal.h"
#include "sx_msort.h"
#include "sx_tpool.h"
#include <stdatomic.h>
#include <sched.h>

double sx_expected_pruned_ltree_nleaves(const sx_stree *st, double b_rate, double d_rate)
{
    return exp((b_rate - d_rate) * (double)sx_stree_gl_height(st, st->root)) * st->n_leaves;
}

typedef struct { double r_bound; int *l_nodes; int n_lnodes; } sx_period;

static int compare_periods(const void *n1, const void *n2)
{
    double v1 = ((const sx_period *)n1)->r_bound, v2 = ((const sx_period *)n2)->r_bound, dif;
    if (v1 - v2 > 0) dif = v1 - v2; else dif = v2 - v1;
    if (dif < DBL_EPSILON) return 0;
    return (v1 > v2) ? 1 : -1;
}

/* ChooseLNodePeriod port: receptor sampling inversely proportional to MRCA distance. */
static int choose_lnode_period(sx_ltree *t, const int *l_pointers, int n_nodes, int t_node, double u_num, int verbosity)
{
    int i = 0, j = 0, done = 0, next_t = 1;
    double *t_times = sx_calloc((size_t)n_nodes, sizeof(double));
    double t_prob = 0, c_prob = 0;
    int w_lnode = -1, w_tnode = t_node;
    int *wl_pointers = sx_malloc(sizeof(int) * (size_t)n_nodes);
    (void)verbosity;
    for (j = 0; j < n_nodes; ++j) t_times[j] = -1;
    memcpy(wl_pointers, l_pointers, sizeof(int) * (size_t)n_nodes);
    while (i < MAX_IT && done == 0) {
        done = 1; next_t = 1;
        for (j = 0; j < n_nodes; ++j) {
            if (t_times[j] != -1) continue;
            w_lnode = wl_pointers[j];
            done = 0;
            if ((t->time[w_lnode] > t->time[w_tnode]) && (w_lnode != w_tnode)) { wl_pointers[j] = t->anc[w_lnode]; next_t = 0; }
            else if (w_lnode == w_tnode) t_times[j] = t->time[w_lnode];
        }
        if (next_t == 1 && done == 0) w_tnode = t->anc[w_tnode];
        ++i;
    }
    for (j = 0; j < n_nodes; ++j) {
        if (t->time[t_node] == t_times[j]) t_times[j] = 0;
        else { t_times[j] = 1 / (t->time[t_node] - t_times[j]); t_prob += t_times[j]; }
    }
    w_lnode = -1;
    for (j = 0; j < n_nodes; ++j) {
        if (u_num > c_prob + (t_times[j] / t_prob)) c_prob += (t_times[j] / t_prob);
        else { w_lnode = l_pointers[j]; break; }
    }
    free(t_times); free(wl_pointers);
    if (w_lnode == -1) sx_check(SX_UNEXPECTED_VALUE, NULL);
    return w_lnode;
}

static void set_from_snode(sx_ltree *lt, int l, const sx_stree *st, int s)
{
    lt->sp_index[l] = st->sp_index[s]; lt->Ne_b[l] = st->Ne_b[s]; lt->mu_mult[l] = st->mu_mult[s]; lt->gtime_mult[l] = st->gtime_mult[s];
}

long sx_sim_bdlh_tree(sx_stree *st, sx_ltree **wlt, int *node_ptrs, double b_rate, double d_rate, double h_rate, double gc_rate, int t_kind,
                      gsl_rng *seed, int min_lleaves, int min_lsleaves, int verbosity, sx_ltree_stats *S, int calc_probs)
{
    sx_ltree *lt = NULL;
    int aux_lnode = -1, anc_lnode = -1, w_lnode = -1, w_lnode2 = -1, w_lnode3 = -1, w_snode = -1;
    int n_leaves = 0, slice_leaves = 0, lt_true_leaves = 0, lt_diffs_true_leaves = 0, diffs_true_leaves = 0, tn_nodes = 0, extra_nodes = 0, n_losses = 0, n_periods = 0, n_ltransf = 0, n_lgc = 0, n_avail_receptors = 0, t_event = 0;
    int next_paralog = 0, node_index = 0, avail_leaves = 0, n_transfer = 0, n_gc = 0, n_nodes = 0;
    double w_prob = d_rate + b_rate + h_rate + gc_rate, dtgc_prob = ((d_rate + h_rate + gc_rate) / (b_rate + d_rate + h_rate + gc_rate)), tgc_prob = ((h_rate + gc_rate) / (b_rate + d_rate + h_rate + gc_rate)), gc_prob = (gc_rate / (b_rate + d_rate + h_rate + gc_rate)), current_ngen = 0, max_ngen = 0, sampled_ngen = 0, rnumber = 0, max_time = 0, gen_time = st->gen_time;
    int bd_iter = 0, l_tree_retries = 0, maxnleaves_reached = 0, done = 0, discarded_t = 0;
    int i, j, k = 0;

    while (done != 1 && S->ntrials <= MAX_IT) {
        if (*wlt != NULL) sx_ltree_free(*wlt);
        *wlt = lt = sx_ltree_new(1, 0, 0, st->max_children, st->gen_time, st->Ne, st->mu);
        sx_stree_reset_siml(st);
        double expected = sx_expected_pruned_ltree_nleaves(st, b_rate, d_rate);
        if (!(expected < (double)MAX_LEAVES)) {
            fprintf(stderr, "\n\nLocus tree birth rate is too hight for the rest of the settings, with an expected number of locus tree leaves bigger than MAX_LEAVES. You could increase this limit changing the environmental variable SIMPHY_MAXLEAVES\n");
            return SX_SETTINGS_ERROR;
        }
        st->l_head[st->root] = lt->root;
        st->n_lnodes[st->root] = 1;
        lt->n_gen[lt->root] = 0; lt->time[lt->root] = 0; lt->conts[lt->root] = st->root;
        next_paralog = 1; n_leaves = 0; lt_true_leaves = 0; lt_diffs_true_leaves = 0; tn_nodes = 1; extra_nodes = 0; n_losses = 0; n_ltransf = 0; n_lgc = 0;
        S->leaves = 0; S->gleaves = 0; S->dups = 0; S->losses = 0; S->transfr = 0; S->gc = 0;
        ++S->ntrials;
        maxnleaves_reached = 0; max_time = 0;

        for (i = 1; i < st->n_nodes; ++i) {
            w_snode = i;
            int anc_s = st->anc[w_snode];
            if (st->n_lnodes[anc_s] == 0) continue;
            int w_lnode3_iter = -2; /* dummy: first iteration takes l_head */
            diffs_true_leaves = 0; slice_leaves = 0;
            for (j = 0; j < st->n_lnodes[anc_s]; ++j) {
                w_lnode3_iter = (j == 0) ? st->l_head[anc_s] : lt->lat[w_lnode3_iter];
                anc_lnode = w_lnode3_iter;
                current_ngen = st->n_gen[anc_s];
                max_ngen = st->n_gen[w_snode];
                sampled_ngen = -log(gsl_rng_uniform_pos(seed)) / (w_prob);
                if (sampled_ngen >= max_ngen - current_ngen) {
                    ++slice_leaves;
                    w_lnode = sx_ltree_add_node(lt);
                    lt->anc[w_lnode] = anc_lnode;
                    SX_CHILD(lt, anc_lnode, lt->n_child[anc_lnode]) = w_lnode;
                    lt->n_child[anc_lnode]++;
                    if (st->n_lnodes[w_snode] == 0) st->l_head[w_snode] = w_lnode;
                    else {
                        w_lnode2 = st->l_head[w_snode];
                        for (k = 1; k < st->n_lnodes[w_snode]; ++k) w_lnode2 = lt->lat[w_lnode2];
                        lt->lat[w_lnode2] = w_lnode;
                    }
                    lt->conts[w_lnode] = w_snode;
                    lt->n_gen[w_lnode] = st->n_gen[w_snode];
                    lt->gen_length[w_lnode] = st->gen_length[w_snode];
                    set_from_snode(lt, w_lnode, st, w_snode);
                    lt->paralog[w_lnode] = lt->paralog[lt->anc[w_lnode]];
                    lt->n_lin[w_lnode] = st->n_replicas[w_snode];
                    lt->n_ilin[w_lnode] = st->n_replicas[w_snode];
                    lt->time[w_lnode] = st->time[w_snode];
                    if (max_time < lt->time[w_lnode]) max_time = lt->time[w_lnode];
                    st->n_lnodes[w_snode]++;
                    ++tn_nodes;
                    if (st->n_child[w_snode] == 0) {
                        ++S->leaves; ++lt_true_leaves; diffs_true_leaves |= 1;
                        S->gleaves += lt->n_lin[w_lnode];
                    }
                } else {
                    w_lnode = sx_ltree_add_node(lt);
                    lt->anc[w_lnode] = anc_lnode;
                    SX_CHILD(lt, anc_lnode, lt->n_child[anc_lnode]) = w_lnode;
                    lt->n_child[anc_lnode]++;
                    lt->conts[w_lnode] = w_snode;
                    lt->paralog[w_lnode] = lt->paralog[lt->anc[w_lnode]];
                    extra_nodes = 0;
                    node_ptrs[0] = w_lnode;
                    n_leaves = 1; avail_leaves = 1; n_nodes = 1;
                    anc_lnode = w_lnode;
                    n_losses = 0; bd_iter = 0;
                    while (current_ngen + sampled_ngen < max_ngen && avail_leaves != 0 && bd_iter <= MAX_IT && slice_leaves < MAX_LEAVES - 1) {
                        ++bd_iter;
                        if (avail_leaves > 1) node_index = (int)gsl_rng_uniform_int(seed, (unsigned long)avail_leaves);
                        else node_index = 0;
                        rnumber = gsl_rng_uniform(seed);
                        if (rnumber > dtgc_prob) {
                            ++slice_leaves;
                            w_lnode = sx_ltree_add_node(lt);
                            w_lnode2 = sx_ltree_add_node(lt);
                            anc_lnode = node_ptrs[node_index];
                            SX_CHILD(lt, anc_lnode, 0) = w_lnode; SX_CHILD(lt, anc_lnode, 1) = w_lnode2;
                            lt->anc[w_lnode] = anc_lnode; lt->anc[w_lnode2] = anc_lnode;
                            lt->n_gen[anc_lnode] = current_ngen + sampled_ngen;
                            lt->gen_length[anc_lnode] = lt->n_gen[anc_lnode] - lt->n_gen[lt->anc[anc_lnode]];
                            lt->time[anc_lnode] = lt->time[lt->anc[anc_lnode]] + lt->gen_length[anc_lnode] * st->gtime_mult[w_snode] * gen_time;
                            lt->kind[anc_lnode] = DUP; lt->n_child[anc_lnode] = 2;
                            set_from_snode(lt, anc_lnode, st, w_snode);
                            lt->conts[anc_lnode] = w_snode;
                            lt->paralog[w_lnode] = lt->paralog[anc_lnode];
                            lt->paralog[w_lnode2] = next_paralog;
                            lt->fmax_nlin[w_lnode2] = 1; lt->n_olin[w_lnode2] = 1;
                            ++next_paralog;
                            node_ptrs[node_index] = w_lnode; node_ptrs[avail_leaves] = w_lnode2;
                            ++avail_leaves; ++n_leaves; n_nodes += 2; ++S->dups;
                        } else if (rnumber > tgc_prob) {
                            --slice_leaves;
                            w_lnode = node_ptrs[node_index];
                            anc_lnode = lt->anc[w_lnode];
                            lt->n_gen[w_lnode] = current_ngen + sampled_ngen;
                            lt->gen_length[w_lnode] = lt->n_gen[w_lnode] - lt->n_gen[anc_lnode];
                            lt->time[w_lnode] = lt->time[lt->anc[w_lnode]] + lt->gen_length[w_lnode] * st->gtime_mult[w_snode] * gen_time;
                            node_ptrs[node_index] = node_ptrs[avail_leaves - 1];
                            --avail_leaves; --n_leaves;
                            lt->sp_index[w_lnode] = st->sp_index[w_snode];
                            lt->kind[w_lnode] = LOSS;
                            lt->Ne_b[w_lnode] = st->Ne_b[w_snode]; lt->mu_mult[w_lnode] = st->mu_mult[w_snode]; lt->gtime_mult[w_lnode] = st->gtime_mult[w_snode];
                            lt->n_lin[w_lnode] = 0; lt->n_ilin[w_lnode] = 0; lt->n_olin[w_lnode] = 0;
                            lt->conts[w_lnode] = w_snode;
                            n_nodes -= 2; extra_nodes += 2; ++n_losses; ++S->losses;
                        } else if (rnumber > gc_prob) {
                            w_lnode = sx_ltree_add_node(lt);
                            anc_lnode = node_ptrs[node_index];
                            SX_CHILD(lt, anc_lnode, 0) = w_lnode;
                            lt->anc[w_lnode] = anc_lnode;
                            lt->n_gen[anc_lnode] = current_ngen + sampled_ngen;
                            lt->gen_length[anc_lnode] = lt->n_gen[anc_lnode] - lt->n_gen[lt->anc[anc_lnode]];
                            lt->time[anc_lnode] = lt->time[lt->anc[anc_lnode]] + lt->gen_length[anc_lnode] * st->gtime_mult[w_snode] * gen_time;
                            lt->kind[anc_lnode] = TRFR; lt->n_child[anc_lnode] = 1;
                            set_from_snode(lt, anc_lnode, st, w_snode);
                            lt->conts[anc_lnode] = w_snode;
                            lt->paralog[w_lnode] = lt->paralog[anc_lnode];
                            node_ptrs[node_index] = w_lnode;
                            ++n_nodes; ++S->transfr;
                        } else {
                            w_lnode = sx_ltree_add_node(lt);
                            anc_lnode = node_ptrs[node_index];
                            SX_CHILD(lt, anc_lnode, 0) = w_lnode;
                            lt->anc[w_lnode] = anc_lnode;
                            lt->n_gen[anc_lnode] = current_ngen + sampled_ngen;
                            lt->gen_length[anc_lnode] = lt->n_gen[anc_lnode] - lt->n_gen[lt->anc[anc_lnode]];
                            lt->time[anc_lnode] = lt->time[lt->anc[anc_lnode]] + lt->gen_length[anc_lnode] * st->gtime_mult[w_snode] * gen_time;
                            lt->kind[anc_lnode] = GC; lt->n_child[anc_lnode] = 1;
                            set_from_snode(lt, anc_lnode, st, w_snode);
                            lt->conts[anc_lnode] = w_snode;
                            lt->paralog[w_lnode] = lt->paralog[anc_lnode];
                            node_ptrs[node_index] = w_lnode;
                            ++n_nodes; ++S->gc;
                        }
                        current_ngen += sampled_ngen;
                        sampled_ngen = -log(gsl_rng_uniform_pos(seed)) / (w_prob * avail_leaves);
                    }
                    if (slice_leaves + 1 >= MAX_LEAVES) { maxnleaves_reached = 1; break; }
                    if (avail_leaves != 0) {
                        if (st->n_lnodes[w_snode] == 0) {
                            w_lnode = node_ptrs[0];
                            lt->n_gen[w_lnode] = st->n_gen[w_snode];
                            lt->gen_length[w_lnode] = lt->n_gen[w_lnode] - lt->n_gen[lt->anc[w_lnode]];
                            lt->time[w_lnode] = st->time[w_snode];
                            set_from_snode(lt, w_lnode, st, w_snode);
                            lt->n_lin[w_lnode] = st->n_replicas[w_snode];
                            lt->n_ilin[w_lnode] = lt->n_lin[w_lnode];
                            lt->conts[w_lnode] = w_snode;
                            if (max_time < lt->time[w_lnode]) max_time = lt->time[w_lnode];
                            st->l_head[w_snode] = w_lnode;
                            aux_lnode = w_lnode;
                            ++st->n_lnodes[w_snode];
                            k = 1;
                        } else {
                            aux_lnode = st->l_head[w_snode];
                            for (k = 1; k < st->n_lnodes[w_snode]; ++k) aux_lnode = lt->lat[aux_lnode];
                            k = 0;
                        }
                    }
                    for (; k < avail_leaves; ++k) {
                        w_lnode = node_ptrs[k];
                        lt->n_gen[w_lnode] = st->n_gen[w_snode];
                        lt->gen_length[w_lnode] = lt->n_gen[w_lnode] - lt->n_gen[lt->anc[w_lnode]];
                        lt->time[w_lnode] = st->time[w_snode];
                        set_from_snode(lt, w_lnode, st, w_snode);
                        lt->n_lin[w_lnode] = st->n_replicas[w_snode];
                        lt->n_ilin[w_lnode] = lt->n_lin[w_lnode];
                        lt->conts[w_lnode] = w_snode;
                        if (max_time < lt->time[w_lnode]) max_time = lt->time[w_lnode];
                        lt->lat[aux_lnode] = w_lnode;
                        aux_lnode = w_lnode;
                        ++st->n_lnodes[w_snode];
                    }
                    if (st->n_child[w_snode] == 0) {
                        S->leaves += n_leaves; lt_true_leaves += n_leaves; diffs_true_leaves |= 1;
                        S->gleaves += n_leaves * st->n_replicas[w_snode];
                    }
                    S->leaves += n_losses;
                    tn_nodes += n_nodes + extra_nodes;
                }
                if (maxnleaves_reached == 1) break;
            }
            if (maxnleaves_reached == 1) break;
            else lt_diffs_true_leaves += diffs_true_leaves;
        }

        if (lt_true_leaves < min_lleaves || maxnleaves_reached == 1 || lt_diffs_true_leaves < min_lsleaves || S->gleaves <= 1) {
            if (verbosity > 2) {
                if (maxnleaves_reached == 1)
                    sx_out("\n\t\tThe locus tree simulation reached the maximum number of locus tree lineages (%d) inside a species tree branch. After checking your simulation parameters carefully (mainly the duplication rate) you can sort this problem out setting the environmental variable SIMPHY_MAXLEAVES with a value higher than %d. Try %d of %d", MAX_LEAVES, MAX_LEAVES, l_tree_retries, MAX_IT);
                else if (lt_true_leaves < min_lleaves)
                    sx_out("\n\t\tLocus tree with %u leaves, less than the minimum %u , restart of the simulation of this tree. Try %d of %d", lt_true_leaves, min_lleaves, l_tree_retries, MAX_IT);
                else if (lt_diffs_true_leaves < min_lsleaves)
                    sx_out("\n\t\tLocus tree with %u leaves from different species, less than the minimum %u , restart of the simulation of this tree. Try %d of %d", lt_diffs_true_leaves, min_lsleaves, l_tree_retries, MAX_IT);
                else
                    sx_out("\n\t\tLocus tree with %u associated gene tree leaves, less than 2, restart of the simulation of this tree. Try %d of %d", S->gleaves, l_tree_retries, MAX_IT);
            }
            ++l_tree_retries;
        } else done = 1;
    }

    if (S->ntrials > MAX_IT) return SX_LOOP_ERROR;

    lt->n_leaves = S->leaves;
    if (lt->n_nodes != tn_nodes) sx_fatal("internal: locus tree node count mismatch (%d vs %d)", lt->n_nodes, tn_nodes);
    lt->n_gleaves = S->gleaves;
    lt->species_tree = st;

    if (S->transfr + S->gc > 0) {
        /* intermediate post-order collapse with room for receptors */
        sx_ltree_collapse(lt, tn_nodes + S->transfr + S->gc, S->leaves, 1);
        lt->n_nodes = tn_nodes;   /* receptors are appended one by one */
        n_periods = tn_nodes - S->leaves + S->losses + 1;
        sx_period *periods = sx_calloc((size_t)n_periods, sizeof(sx_period));
        periods[0].r_bound = max_time; periods[0].l_nodes = NULL; periods[0].n_lnodes = 0;
        j = 1;
        int *avail_receptors = sx_calloc((size_t)S->leaves, sizeof(int));
        for (i = 0; i < tn_nodes; ++i) {
            if (lt->kind[i] == SP && lt->n_child[i] == 0) continue;
            if (j >= n_periods) sx_fatal("internal: period overflow");
            sx_period *wp = &periods[j];
            wp->l_nodes = sx_calloc((size_t)S->leaves, sizeof(int));
            wp->n_lnodes = lt->n_child[i];
            for (k = 0; k < lt->n_child[i]; ++k) wp->l_nodes[k] = SX_CHILD(lt, i, k);
            wp->r_bound = lt->time[i];
            ++j;
        }
        sx_msort(periods, (size_t)n_periods, sizeof(sx_period), compare_periods);
        for (i = 0; i < n_periods; ++i) {
            sx_period *wp = &periods[i];
            sx_period *wp2 = (i + 1 < n_periods) ? &periods[i + 1] : NULL;
            double next_bound = wp2 ? wp2->r_bound : INFINITY;
            n_avail_receptors = 0; discarded_t = 0; w_lnode2 = -1;
            for (j = 0; j < wp->n_lnodes; ++j) {
                w_lnode = wp->l_nodes[j];
                if (lt->time[w_lnode] > next_bound) {
                    if (!wp2) sx_fatal("internal: period sweep overflow");
                    wp2->l_nodes[wp2->n_lnodes] = w_lnode; wp2->n_lnodes += 1;
                } else if (wp->n_lnodes > 1 && ((lt->kind[w_lnode] == TRFR || lt->kind[w_lnode] == GC) && lt->time[w_lnode] == next_bound)) {
                    w_snode = lt->conts[w_lnode]; w_lnode2 = w_lnode; t_event = lt->kind[w_lnode];
                } else if (wp->n_lnodes == 1 && ((lt->kind[w_lnode] == TRFR || lt->kind[w_lnode] == GC) && lt->time[w_lnode] == next_bound)) {
                    w_lnode2 = w_lnode; discarded_t = 1;
                }
            }
            if (w_lnode2 != -1) {
                int did_transfer = 0;
                if (discarded_t == 0) {
                    for (j = 0; j < wp->n_lnodes; ++j) {
                        w_lnode = wp->l_nodes[j];
                        switch (t_event) {
                        case TRFR:
                            if (lt->conts[w_lnode] != w_snode) avail_receptors[n_avail_receptors++] = w_lnode;
                            break;
                        case GC:
                            if (lt->conts[w_lnode] == w_snode && w_lnode != w_lnode2) avail_receptors[n_avail_receptors++] = w_lnode;
                            break;
                        default: sx_check(SX_UNEXPECTED_VALUE, NULL); break;
                        }
                    }
                    if (n_avail_receptors > 0) {
                        switch (t_kind) {
                        case 1: w_lnode = w_lnode2; w_lnode2 = choose_lnode_period(lt, avail_receptors, n_avail_receptors, w_lnode2, gsl_rng_uniform_pos(seed), verbosity); break;
                        case 0: w_lnode = w_lnode2; w_lnode2 = avail_receptors[gsl_rng_uniform_int(seed, (unsigned long)n_avail_receptors)]; break;
                        default: sx_check(SX_SETTINGS_ERROR, NULL); break;
                        }
                        w_lnode3 = sx_ltree_add_node(lt);
                        if (w_lnode3 != tn_nodes + n_transfer + n_gc) sx_fatal("internal: receptor id mismatch");
                        anc_lnode = lt->anc[w_lnode2];
                        lt->anc[w_lnode3] = anc_lnode;
                        lt->paralog[w_lnode3] = lt->paralog[w_lnode2];
                        lt->conts[w_lnode3] = lt->conts[w_lnode2];
                        lt->mu_mult[w_lnode3] = lt->mu_mult[w_lnode2];
                        lt->gtime_mult[w_lnode3] = lt->gtime_mult[w_lnode2];
                        lt->Ne_b[w_lnode3] = lt->Ne_b[w_lnode2];
                        lt->time[w_lnode3] = lt->time[w_lnode];
                        lt->n_gen[w_lnode3] = lt->n_gen[anc_lnode] + (lt->time[w_lnode3] - lt->time[anc_lnode]) / lt->gtime_mult[w_lnode3] / gen_time;
                        lt->gen_length[w_lnode3] = lt->n_gen[w_lnode3] - lt->n_gen[anc_lnode];
                        lt->sp_index[w_lnode3] = st->sp_index[lt->conts[w_lnode3]];
                        for (k = 0; k < lt->n_child[anc_lnode]; ++k)
                            if (SX_CHILD(lt, anc_lnode, k) == w_lnode2) SX_CHILD(lt, anc_lnode, k) = w_lnode3;
                        lt->gen_length[w_lnode2] = lt->n_gen[w_lnode2] - lt->n_gen[w_lnode3];
                        lt->anc[w_lnode2] = w_lnode;
                        lt->fmax_nlin[w_lnode2] = 1; lt->n_olin[w_lnode2] = 1;
                        SX_CHILD(lt, w_lnode, 1) = w_lnode2;
                        lt->n_child[w_lnode] = 2;
                        if (t_event == TRFR) { lt->kind[w_lnode3] = RTRFR; ++n_transfer; }
                        else { lt->kind[w_lnode3] = RGC; ++n_gc; }
                        did_transfer = 1;
                    }
                }
                if (!did_transfer) {
                    /* discarded: bypass the donor node */
                    int child0 = SX_CHILD(lt, w_lnode2, 0), anc2 = lt->anc[w_lnode2];
                    lt->anc[child0] = anc2;
                    for (k = 0; k < lt->n_child[anc2]; ++k)
                        if (SX_CHILD(lt, anc2, k) == w_lnode2) { SX_CHILD(lt, anc2, k) = child0; break; }
                    lt->gen_length[child0] = lt->n_gen[child0] - lt->n_gen[anc2];
                    if (lt->kind[w_lnode2] == TRFR) ++n_ltransf; else ++n_lgc;
                }
            }
        }
        for (i = 0; i < n_periods; ++i) free(periods[i].l_nodes);
        free(periods); free(avail_receptors);
        lt->root = tn_nodes - 1;
    }
    S->transfr -= n_ltransf;
    S->gc -= n_lgc;
    st->locus_tree = lt;

    if (S->dups > 0 || S->transfr > 0 || S->gc > 0) {
        sx_ltree_collapse(lt, tn_nodes + S->transfr - n_ltransf + S->gc - n_lgc, S->leaves + S->transfr + S->gc, 1);
        if (calc_probs) sx_calc_probs_nlineages(lt, lt->root, lt->Ne, 0, verbosity);
    } else
        sx_ltree_collapse(lt, tn_nodes + S->transfr - n_ltransf + S->gc - n_lgc, S->leaves + S->transfr + S->gc, 1);
    return SX_NO_ERROR;
}

/* ---------------- CalcProbsNLineagesLTree (compat) ---------------- */
static sx_tpool *calc_pool = NULL;
void sx_calc_probs_set_pool(sx_tpool *pool) { calc_pool = pool; }

typedef struct { double *o, *ip; int fmax, i0, i1; double t; int Ne; atomic_int *left; } row_task;
static void row_task_run(void *vp)
{
    row_task *rt = vp;
    int stride = rt->i1 < 0 ? -rt->i1 : 1, end = rt->i1 < 0 ? rt->fmax : rt->i1 - 1;
    for (int i = rt->i0; i <= end; i += stride) {
        double comp = 0;
        for (int j = i; j <= rt->fmax; ++j)
            sx_skahan_sum(&rt->o[i], sx_kahan_logscale_prob_coal(j, i, rt->t, rt->Ne) * rt->ip[j], &comp);
    }
    atomic_fetch_sub(rt->left, 1);
    free(rt);
}
static void calc_rows_parallel(double *o, double *ip, int fmax, double t, int Ne)
{
    /* rows cost ~ (fmax-i): interleave so chunks are balanced; chunk = 1 row for large i, more for small */
    atomic_int left = 0;
    int nthr = sx_tpool_size(calc_pool) + 1;
    int chunks = 4 * nthr; if (chunks > fmax) chunks = fmax;
    /* assign rows round-robin to chunks: chunk c gets rows i with (i-1) % chunks == c → balanced */
    for (int c = 0; c < chunks; ++c) {
        /* emulate round robin with contiguous tasks: rows i = c+1, c+1+chunks, ... handled by a loop inside one task */
        row_task *rt = sx_malloc(sizeof(*rt));
        rt->o = o; rt->ip = ip; rt->fmax = fmax; rt->t = t; rt->Ne = Ne; rt->left = &left;
        rt->i0 = c + 1; rt->i1 = -chunks;   /* stride encoded: negative i1 means stride */
        atomic_fetch_add(&left, 1);
        sx_tpool_submit_front(calc_pool, row_task_run, rt);
    }
    while (atomic_load(&left) > 0) { if (!sx_tpool_help(calc_pool)) sched_yield(); }
}
/* Debug: log every lineage-count probability evaluation (i, o, t, Ne, value) to $SIMPHY_X_PROBLOG. */
static FILE *problog = NULL;
static int problog_checked = 0;
static double kahan_logged(int i, int o, double t, int Ne, const char *where)
{
    double p = sx_kahan_logscale_prob_coal(i, o, t, Ne);
    if (!problog_checked) { problog_checked = 1; const char *f = getenv("SIMPHY_X_PROBLOG"); if (f) problog = fopen(f, "w"); }
    if (problog) fprintf(problog, "%s %d %d %.17g %d %.17g\n", where, i, o, t, Ne, p);
    return p;
}
int sx_calc_probs_nlineages(sx_ltree *t, int node, int Ne, int in_subtree, int verbosity)
{
    int *max_nlin = NULL, i = 0, j = 0, k = 0, fmax_nlin = 0, combfmax_nlin = 1, return_nlin = 0, nbranch_wlin = 0;
    double comp = 0, comp_opsum = 0, o_prob_sum = 0, bkp_prob = 0;
    int nch = t->n_child[node];
    int kind = t->kind[node];

    if (nch != 0) {
        max_nlin = sx_calloc((size_t)nch, sizeof(int));
        for (i = 0; i < nch; ++i) {
            int c = SX_CHILD(t, node, i);
            if (i == 1 && (kind == DUP || kind == TRFR || kind == GC)) max_nlin[i] = sx_calc_probs_nlineages(t, c, Ne, 1, verbosity);
            else max_nlin[i] = sx_calc_probs_nlineages(t, c, Ne, in_subtree, verbosity);
            fmax_nlin += max_nlin[i];
            combfmax_nlin *= max_nlin[i];
            if (max_nlin[i] != 0) ++nbranch_wlin;
        }
        if (in_subtree == 1) {
            free(t->i_probs[node]); free(t->i_combprobs[node]); free(t->o_probs[node]);
            t->i_probs[node] = sx_calloc((size_t)fmax_nlin + 1, sizeof(double));
            t->i_combprobs[node] = sx_calloc((size_t)(combfmax_nlin > 0 ? combfmax_nlin : 1), sizeof(double));
            t->o_probs[node] = sx_calloc((size_t)fmax_nlin + 1, sizeof(double));
            int c0 = SX_CHILD(t, node, 0), c1 = nch > 1 ? SX_CHILD(t, node, 1) : -1;
            if (nch > 1 && max_nlin[1] == 1) { bkp_prob = t->o_probs[c1][1]; t->o_probs[c1][1] = 1; }
            switch (nbranch_wlin) {
            case 2:
                k = 0;
                for (i = 2; i <= fmax_nlin; ++i)
                    for (j = 1; j < i; ++j) {
                        if (j <= max_nlin[0] && (i - j) <= max_nlin[1]) {
                            t->i_combprobs[node][k] = t->o_probs[c0][j] * t->o_probs[c1][i - j];
                            t->i_probs[node][i] += t->i_combprobs[node][k];
                            ++k;
                        }
                    }
                break;
            case 1:
                for (j = 0; j < nch; ++j) {
                    if (max_nlin[j] != 0) {
                        int cj = SX_CHILD(t, node, j);
                        for (i = 1; i <= fmax_nlin; ++i) t->i_probs[node][i] = t->o_probs[cj][i];
                        break;
                    }
                }
                break;
            case 0:
                fmax_nlin = 0;
                t->i_probs[node][0] = 1;
                t->o_probs[node][0] = 1;
                break;
            default:
                sx_check(SX_UNEXPECTED_VALUE, "Not implemented yet\n");
                break;
            }
        }
    } else {
        fmax_nlin = t->n_ilin[node];
        free(t->i_probs[node]); free(t->o_probs[node]);
        t->i_probs[node] = sx_calloc((size_t)fmax_nlin + 1, sizeof(double));
        t->o_probs[node] = sx_calloc((size_t)fmax_nlin + 1, sizeof(double));
        t->i_probs[node][fmax_nlin] = 1;
    }
    return_nlin = fmax_nlin;

    if (in_subtree == 1) {
        int anc = t->anc[node];
        if (anc == -1 || fmax_nlin == 0) {
            t->n_olin[node] = 0;
            t->o_probs[node][0] = 1;
            return_nlin = 0;
        } else if (fmax_nlin == 1) {
            t->n_olin[node] = 1;
            t->o_probs[node][1] = 1;
        } else {
            comp_opsum = 0;
            int NeN = t->Ne_b[node] == 0 ? Ne : t->Ne_b[node];
            if (calc_pool && fmax_nlin >= 24 && !problog) {
                /* parallel rows, then the sequential break rule (rows after the break are zeroed as SimPhy leaves them) */
                calc_rows_parallel(t->o_probs[node], t->i_probs[node], fmax_nlin, t->gen_length[node], NeN);
                for (i = 1; i <= fmax_nlin; ++i) {
                    sx_skahan_sum(&o_prob_sum, t->o_probs[node][i], &comp_opsum);
                    if (o_prob_sum >= 1) break;
                }
                for (j = i + 1; j <= fmax_nlin; ++j) t->o_probs[node][j] = 0;
            } else {
                for (i = 1; i <= fmax_nlin; ++i) {
                    comp = 0;
                    for (j = i; j <= fmax_nlin; ++j)
                        sx_skahan_sum(&t->o_probs[node][i], kahan_logged(j, i, t->gen_length[node], NeN, "C") * t->i_probs[node][j], &comp);
                    sx_skahan_sum(&o_prob_sum, t->o_probs[node][i], &comp_opsum);
                    if (o_prob_sum >= 1) break;
                }
            }
            if (anc != -1 && ((t->kind[anc] == DUP || t->kind[anc] == TRFR || t->kind[anc] == GC) && SX_CHILD(t, anc, 1) == node)) {
                t->n_olin[node] = 1;
                return_nlin = 1;
            }
        }
        if (nch > 0 && max_nlin[1] == 1) t->o_probs[SX_CHILD(t, node, 1)][1] = bkp_prob;
    }
    free(max_nlin);
    t->fmax_nlin[node] = return_nlin;
    return return_nlin;
}

/* ---------------- SampleNLineagesCDFLTree (compat) ---------------- */
void sx_sample_nlineages_cdf(sx_ltree *t, int node, int n_gleaves, int Ne, int verbosity, gsl_rng *seed)
{
    int i = 0, j = 0, prob_padding = 0, prob_combs = 0, *mapcomb_lin0 = NULL, n_branch_wlin = 0, wn_olin0 = 0, wn_olin1 = 0, child = 0;
    double p = 0, u = 0;
    int kt = 0, ktmax = 0;
    double sum = 0, cden = 0, comp = 0;
    int nch = t->n_child[node];
    if (nch == 0) return;
    if (t->n_olin[node] == 0) {
        for (i = 0; i < nch; ++i) sx_sample_nlineages_cdf(t, SX_CHILD(t, node, i), n_gleaves, Ne, verbosity, seed);
        return;
    }
    int c0 = SX_CHILD(t, node, 0), c1 = SX_CHILD(t, node, 1);
    u = gsl_rng_uniform(seed);
    ktmax = t->fmax_nlin[c0] + t->fmax_nlin[c1];
    cden = t->o_probs[node][t->n_olin[node]];
    comp = 0;
    kt = t->n_olin[node] - 1;
    if (t->n_olin[node] < ktmax) {
        do {
            ++kt;
            sx_skahan_sum(&sum, kahan_logged(kt, t->n_olin[node], t->gen_length[node], t->Ne_b[node] != 0 ? t->Ne_b[node] : Ne, "S") * t->i_probs[node][kt], &comp);
            p = sum / cden;
        } while (p < u && kt < ktmax);
    } else if (t->n_olin[node] == ktmax) kt = ktmax;
    else sx_check(SX_UNEXPECTED_VALUE, "");

    for (i = 0; i < nch; ++i) {
        int c = SX_CHILD(t, node, i);
        if (t->o_probs[c][0] == 0) { ++n_branch_wlin; child = i; }
    }
    switch (n_branch_wlin) {
    case 2:
        t->n_ilin[node] = kt;
        prob_padding = 0; prob_combs = 0;
        mapcomb_lin0 = sx_calloc((size_t)kt, sizeof(int));
        for (i = 2; i < kt; ++i)
            for (j = 1; j < i; ++j)
                if (j <= t->fmax_nlin[c0] && (i - j) <= t->fmax_nlin[c1]) ++prob_padding;
        for (j = 1; j < kt; ++j)
            if (j <= t->fmax_nlin[c0] && (kt - j) <= t->fmax_nlin[c1]) { mapcomb_lin0[prob_combs] = j; ++prob_combs; }
        wn_olin0 = mapcomb_lin0[sx_sample_ndoubles((size_t)prob_combs, t->i_combprobs[node] + prob_padding, seed)];
        wn_olin1 = kt - wn_olin0;
        t->n_olin[c0] = wn_olin0;
        t->n_olin[c1] = wn_olin1;
        free(mapcomb_lin0);
        sx_sample_nlineages_cdf(t, c0, n_gleaves, Ne, verbosity, seed);
        sx_sample_nlineages_cdf(t, c1, n_gleaves, Ne, verbosity, seed);
        break;
    case 1:
        t->n_ilin[node] = kt;
        t->n_olin[SX_CHILD(t, node, child)] = kt;
        sx_sample_nlineages_cdf(t, SX_CHILD(t, node, child), n_gleaves, Ne, verbosity, seed);
        break;
    case 0:
        break;
    default:
        sx_check(SX_UNEXPECTED_VALUE, ": not implemented yet\n");
        break;
    }
}

/* ======================= fast mode ======================= */
#include "sx_unif.h"

int sx_calc_probs_nlineages_fast(sx_ltree *t, int node, int Ne, int in_subtree, sx_unif_ws *ws)
{
    int *max_nlin = NULL, i, j, k, fmax_nlin = 0, combfmax_nlin = 1, return_nlin = 0, nbranch_wlin = 0;
    double bkp_prob = 0;
    int nch = t->n_child[node], kind = t->kind[node];
    if (nch != 0) {
        max_nlin = sx_calloc((size_t)nch, sizeof(int));
        for (i = 0; i < nch; ++i) {
            int c = SX_CHILD(t, node, i);
            if (i == 1 && (kind == DUP || kind == TRFR || kind == GC)) max_nlin[i] = sx_calc_probs_nlineages_fast(t, c, Ne, 1, ws);
            else max_nlin[i] = sx_calc_probs_nlineages_fast(t, c, Ne, in_subtree, ws);
            fmax_nlin += max_nlin[i];
            combfmax_nlin *= max_nlin[i];
            if (max_nlin[i] != 0) ++nbranch_wlin;
        }
        if (in_subtree == 1) {
            free(t->i_probs[node]); free(t->i_combprobs[node]); free(t->o_probs[node]);
            t->i_probs[node] = sx_calloc((size_t)fmax_nlin + 1, sizeof(double));
            t->i_combprobs[node] = NULL;   /* split probabilities are recomputed at sampling time */
            t->o_probs[node] = sx_calloc((size_t)fmax_nlin + 1, sizeof(double));
            int c0 = SX_CHILD(t, node, 0), c1 = nch > 1 ? SX_CHILD(t, node, 1) : -1;
            if (nch > 1 && max_nlin[1] == 1) { bkp_prob = t->o_probs[c1][1]; t->o_probs[c1][1] = 1; }
            switch (nbranch_wlin) {
            case 2:
                k = 0;
                for (i = 2; i <= fmax_nlin; ++i)
                    for (j = 1; j < i; ++j)
                        if (j <= max_nlin[0] && (i - j) <= max_nlin[1]) {
                            t->i_probs[node][i] += t->o_probs[c0][j] * t->o_probs[c1][i - j];
                            ++k;
                        }
                break;
            case 1:
                for (j = 0; j < nch; ++j)
                    if (max_nlin[j] != 0) {
                        int cj = SX_CHILD(t, node, j);
                        for (i = 1; i <= fmax_nlin; ++i) t->i_probs[node][i] = t->o_probs[cj][i];
                        break;
                    }
                break;
            case 0:
                fmax_nlin = 0; t->i_probs[node][0] = 1; t->o_probs[node][0] = 1;
                break;
            default:
                sx_check(SX_UNEXPECTED_VALUE, "Not implemented yet\n");
            }
        }
    } else {
        fmax_nlin = t->n_ilin[node];
        free(t->i_probs[node]); free(t->o_probs[node]);
        t->i_probs[node] = sx_calloc((size_t)fmax_nlin + 1, sizeof(double));
        t->o_probs[node] = sx_calloc((size_t)fmax_nlin + 1, sizeof(double));
        t->i_probs[node][fmax_nlin] = 1;
    }
    return_nlin = fmax_nlin;
    if (in_subtree == 1) {
        int anc = t->anc[node];
        if (anc == -1 || fmax_nlin == 0) { t->n_olin[node] = 0; t->o_probs[node][0] = 1; return_nlin = 0; }
        else if (fmax_nlin == 1) { t->n_olin[node] = 1; t->o_probs[node][1] = 1; }
        else {
            /* o_probs = i_probs transported through the branch */
            sx_unif_forward(ws, fmax_nlin, t->gen_length[node], (double)(t->Ne_b[node] == 0 ? Ne : t->Ne_b[node]), t->i_probs[node], t->o_probs[node]);
            if (anc != -1 && ((t->kind[anc] == DUP || t->kind[anc] == TRFR || t->kind[anc] == GC) && SX_CHILD(t, anc, 1) == node)) { t->n_olin[node] = 1; return_nlin = 1; }
        }
        if (nch > 0 && max_nlin[1] == 1) t->o_probs[SX_CHILD(t, node, 1)][1] = bkp_prob;
    }
    free(max_nlin);
    t->fmax_nlin[node] = return_nlin;
    return return_nlin;
}

void sx_sample_nlineages_fast(sx_ltree *t, int node, int Ne, gsl_rng *seed, sx_unif_ws *ws)
{
    int nch = t->n_child[node];
    if (nch == 0) return;
    if (t->n_olin[node] == 0) {
        for (int i = 0; i < nch; ++i) sx_sample_nlineages_fast(t, SX_CHILD(t, node, i), Ne, seed, ws);
        return;
    }
    int c0 = SX_CHILD(t, node, 0), c1 = SX_CHILD(t, node, 1);
    int o = t->n_olin[node];
    int ktmax = t->fmax_nlin[c0] + t->fmax_nlin[c1];
    int kt;
    double u = gsl_rng_uniform(seed);
    if (o < ktmax) {
        /* P(kt | o) ∝ P_{kt,o}(t) i_probs[kt] for kt in [o..ktmax]; column by uniformisation */
        int n = ktmax - o + 1;
        double *col = sx_malloc(sizeof(double) * (size_t)n);
        double Nn = (double)(t->Ne_b[node] != 0 ? t->Ne_b[node] : Ne);
        double LT = (double)ktmax * (double)(ktmax - 1) / (2.0 * Nn) * t->gen_length[node];
        double tot = 0;
        if (LT <= SX_UNIF_MAX_LAMBDA_T) {
            sx_unif_column(ws, ktmax, o, t->gen_length[node], Nn, col);
            for (int k = 0; k < n; ++k) { col[k] *= t->i_probs[node][o + k]; tot += col[k]; }
        } else {
            /* many expected events: posterior weights decay super-exponentially in kt;
             * evaluate the Tavaré series in increasing kt and stop once negligible */
            double mx = 0; int small_run = 0;
            for (int k = 0; k < n; ++k) {
                double ip = t->i_probs[node][o + k];
                double neglig = (mx > 0 && ip > 0) ? 1e-22 * mx / ip : 0.0;
                double w = ip > 0 ? sx_tavare_prob_bounded(ws, o + k, o, t->gen_length[node], Nn, neglig) * ip : 0.0;
                col[k] = w; tot += w;
                if (w > mx) mx = w;
                if (mx > 0 && w < 1e-18 * mx) { if (++small_run >= 4) { for (int kk = k + 1; kk < n; ++kk) col[kk] = 0; break; } }
                else small_run = 0;
            }
        }
        kt = ktmax;
        if (tot > 0) {
            double acc = 0, target = u * tot;
            for (int k = 0; k < n; ++k) { acc += col[k]; if (acc >= target) { kt = o + k; break; } }
        }
        free(col);
    } else if (o == ktmax) kt = ktmax;
    else { sx_check(SX_UNEXPECTED_VALUE, ""); return; }

    int n_branch_wlin = 0, child = 0;
    for (int i = 0; i < nch; ++i) { int c = SX_CHILD(t, node, i); if (t->o_probs[c][0] == 0) { ++n_branch_wlin; child = i; } }
    switch (n_branch_wlin) {
    case 2: {
        t->n_ilin[node] = kt;
        int prob_combs = 0;
        int *map0 = sx_calloc((size_t)kt, sizeof(int));
        for (int j = 1; j < kt; ++j)
            if (j <= t->fmax_nlin[c0] && (kt - j) <= t->fmax_nlin[c1]) { map0[prob_combs] = j; ++prob_combs; }
        /* the split weights are the products o_probs[c0][j] * o_probs[c1][kt-j] (same values and order as i_combprobs) */
        double *cw = sx_malloc(sizeof(double) * (size_t)(prob_combs > 0 ? prob_combs : 1));
        for (int q = 0; q < prob_combs; ++q) cw[q] = t->o_probs[c0][map0[q]] * t->o_probs[c1][kt - map0[q]];
        int wn0 = map0[sx_sample_ndoubles((size_t)prob_combs, cw, seed)];
        free(cw);
        t->n_olin[c0] = wn0; t->n_olin[c1] = kt - wn0;
        free(map0);
        sx_sample_nlineages_fast(t, c0, Ne, seed, ws);
        sx_sample_nlineages_fast(t, c1, Ne, seed, ws);
        break;
    }
    case 1:
        t->n_ilin[node] = kt;
        t->n_olin[SX_CHILD(t, node, child)] = kt;
        sx_sample_nlineages_fast(t, SX_CHILD(t, node, child), Ne, seed, ws);
        break;
    case 0: break;
    default: sx_check(SX_UNEXPECTED_VALUE, ": not implemented yet\n");
    }
}
