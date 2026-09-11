#include "sx_stree.h"
#include "sx_msort.h"

static void add_constant_ngen(sx_stree *t, int node, double constant, double gen_time)
{
    t->n_gen[node] += constant;
    t->time[node] = t->time[t->anc[node]] + t->gen_length[node] * t->gtime_mult[node] * gen_time;
    for (int i = 0; i < t->n_child[node]; ++i) add_constant_ngen(t, SX_CHILD(t, node, i), constant, gen_time);
}

long sx_new_bds_tree(sx_stree **out_tree, int leaves, double time, double b_rate, double d_rate, double gen_time, int Ne, double mu,
                     int ind_per_sp, double outgroup, int complete, int mrca_time, gsl_rng *seed, int verbosity)
{
    int eq_rates = 0, j = 0, k = 0, node_index = 0, n_leaves = 0, pn_leaves = 0, pn_events = 0, n_nodes = 0, avail_leaves = 0, extra_nodes = 0, iter = 0;
    double w_prob = d_rate + b_rate, b_prob = (b_rate / (d_rate + b_rate)), random = 0, inv_leaves = 0, res1 = 0, res2 = 0, *i_nodes = NULL, s_time = 0, c_time = 0, r_time = 0;
    int w_node1 = -1, w_node2 = -1, anc_node = -1, *node_ptrs = NULL;
    sx_stree *t = NULL;

    if (w_prob <= 0 || (b_rate < d_rate && leaves != 0)) return SX_SETTINGS_ERROR;
    if (gen_time == 0) return SX_UNEXPECTED_VALUE;
    if (d_rate == b_rate) eq_rates = 1;

    if (leaves > 0) {
        /* ---- BDSA ---- */
        if (verbosity > 3) printf("\n\tUsing BDSA algorithm (fixed number of leaves) to simulate the trees...");
        if (outgroup > 0) t = sx_stree_new(((leaves + 1) * 2) - 1, leaves + 1, leaves + 1, 2, gen_time, Ne, mu);
        else t = sx_stree_new((leaves * 2) - 1, leaves, leaves, 2, gen_time, Ne, mu);
        *out_tree = t;

        i_nodes = sx_calloc((size_t)(leaves - 1), sizeof(double));
        node_ptrs = sx_calloc((size_t)t->n_nodes, sizeof(int));

        if (time == 0) {
            random = gsl_rng_uniform_pos(seed);
            inv_leaves = 1 / (double)leaves;
            if (eq_rates > 0) time = 1 / (b_rate * (pow(random, -inv_leaves) - 1));
            else {
                res1 = pow(random, inv_leaves);
                time = log((1 - (res1 * d_rate / b_rate)) / (1 - res1)) / (b_rate - d_rate);
            }
            for (j = 0; j < leaves - 1; ++j) {
                random = gsl_rng_uniform_pos(seed);
                if (eq_rates > 0) i_nodes[j] = (random * time) / (1 + (b_rate * time * (1 - random)));
                else {
                    res1 = exp((b_rate - d_rate) * -time);
                    res2 = b_rate - (d_rate * res1);
                    i_nodes[j] = (log((res2 - (d_rate * random * (1 - res1))) / (res2 - (b_rate * random * (1 - res1)))) / (b_rate - d_rate));
                }
            }
            sx_msort(i_nodes, (size_t)(leaves - 1), sizeof(double), sx_compare_dbl);
            r_time = i_nodes[leaves - 2];
            for (j = 0; j < t->n_nodes; ++j) {
                node_ptrs[j] = j;
                t->n_gen[j] = (r_time + r_time * outgroup) / gen_time;
                t->time[j] = r_time + r_time * outgroup;
            }
            n_leaves = leaves;
            for (j = 0; j < leaves - 1; ++j) {
                node_index = (int)gsl_rng_uniform_int(seed, (unsigned long)n_leaves);
                w_node1 = node_ptrs[node_index];
                --n_leaves;
                node_ptrs[node_index] = node_ptrs[n_leaves];
                if (n_leaves > 1) node_index = (int)gsl_rng_uniform_int(seed, (unsigned long)n_leaves);
                else node_index = 0;
                w_node2 = node_ptrs[node_index];
                node_ptrs[node_index] = j + leaves;
                anc_node = j + leaves;
                SX_CHILD(t, anc_node, 0) = w_node1; SX_CHILD(t, anc_node, 1) = w_node2;
                t->anc[w_node1] = anc_node; t->anc[w_node2] = anc_node;
                t->time[anc_node] = r_time + r_time * outgroup - i_nodes[j];
                t->n_gen[anc_node] = t->time[anc_node] / gen_time;
                t->n_child[anc_node] = 2;
                t->gen_length[w_node1] = t->n_gen[w_node1] - t->n_gen[anc_node];
                t->gen_length[w_node2] = t->n_gen[w_node2] - t->n_gen[anc_node];
            }
            t->root = anc_node;
        } else {
            /* TMRCA given */
            r_time = time;
            n_leaves = (int)gsl_rng_uniform_int(seed, (unsigned long)(leaves - 1)) + 1;
            for (j = 0; j < t->n_nodes; ++j) {
                node_ptrs[j] = j;
                t->n_gen[j] = (time + time * outgroup) / gen_time;
                t->time[j] = time + time * outgroup;
            }
            if (leaves == 2) k = 3;
            else if (n_leaves > 1 && n_leaves < leaves - 1) { k = 1; pn_leaves = 0; }
            else if (n_leaves == 1) { pn_leaves = 1; n_leaves = leaves - pn_leaves; k = 2; }
            else { k = 1; pn_leaves = 0; }
            pn_events = 0;
            for (; k <= 2; ++k) {
                extra_nodes = pn_leaves == 0 ? 0 : pn_leaves * 2 - 1;
                (void)extra_nodes;
                pn_events = pn_leaves == 0 ? 0 : pn_leaves - 1;
                for (j = 0; j < n_leaves - 1; ++j) {
                    random = gsl_rng_uniform_pos(seed);
                    if (eq_rates > 0) i_nodes[pn_events + j] = (random * time) / (1 + (b_rate * time * (1 - random)));
                    else {
                        res1 = exp((b_rate - d_rate) * -time);
                        res2 = b_rate - (d_rate * res1);
                        i_nodes[pn_events + j] = (log((res2 - (d_rate * random * (1 - res1))) / (res2 - (b_rate * random * (1 - res1)))) / (b_rate - d_rate));
                    }
                }
                sx_msort(i_nodes + pn_events, (size_t)(n_leaves - 1), sizeof(double), sx_compare_dbl);
                avail_leaves = n_leaves;
                for (j = 0; j < n_leaves - 1; ++j) {
                    node_index = (int)gsl_rng_uniform_int(seed, (unsigned long)avail_leaves) + pn_leaves;
                    w_node1 = node_ptrs[node_index];
                    --avail_leaves;
                    node_ptrs[node_index] = node_ptrs[avail_leaves + pn_leaves];
                    if (avail_leaves > 1) node_index = (int)gsl_rng_uniform_int(seed, (unsigned long)avail_leaves) + pn_leaves;
                    else node_index = 0 + pn_leaves;
                    w_node2 = node_ptrs[node_index];
                    node_ptrs[node_index] = j + leaves + pn_events;
                    anc_node = j + leaves + pn_events;
                    SX_CHILD(t, anc_node, 0) = w_node1; SX_CHILD(t, anc_node, 1) = w_node2;
                    t->anc[w_node1] = anc_node; t->anc[w_node2] = anc_node;
                    t->time[anc_node] = r_time + r_time * outgroup - i_nodes[j + pn_events];
                    t->n_gen[anc_node] = t->time[anc_node] / gen_time;
                    t->n_child[anc_node] = 2;
                    t->gen_length[w_node1] = t->n_gen[w_node1] - t->n_gen[anc_node];
                    t->gen_length[w_node2] = t->n_gen[w_node2] - t->n_gen[anc_node];
                }
                pn_leaves = n_leaves;
                n_leaves = leaves - pn_leaves;
            }
            int root = leaves * 2 - 2;
            t->root = root;
            t->time[root] = r_time * outgroup;
            t->gen_length[root] = 0;
            t->n_gen[root] = t->time[root] / gen_time;
            t->n_child[root] = 2;
            if (pn_leaves == 1) { SX_CHILD(t, root, 0) = leaves * 2 - 3; SX_CHILD(t, root, 1) = leaves - 1; }
            else if (n_leaves == 1) { SX_CHILD(t, root, 0) = 0; SX_CHILD(t, root, 1) = leaves * 2 - 3; }
            else { SX_CHILD(t, root, 0) = leaves + n_leaves - 2; SX_CHILD(t, root, 1) = leaves * 2 - 3; }
            int c1 = SX_CHILD(t, root, 1), c0 = SX_CHILD(t, root, 0);
            t->anc[c1] = root; t->gen_length[c1] = t->n_gen[c1] - t->n_gen[root];
            t->anc[c0] = root; t->gen_length[c0] = t->n_gen[c0] - t->n_gen[root];
        }

        for (j = 0; j < leaves; ++j) { t->sp_index[j] = j + 1; t->n_replicas[j] = ind_per_sp; }
        if (outgroup > 0) {
            w_node1 = leaves * 2;      /* outgroup leaf */
            w_node2 = leaves * 2 - 1;  /* new root */
            t->anc[w_node1] = w_node2;
            SX_CHILD(t, w_node2, 0) = t->root; SX_CHILD(t, w_node2, 1) = w_node1;
            t->n_gen[w_node2] = 0; t->time[w_node2] = 0; t->n_child[w_node2] = 2;
            t->gen_length[t->root] = outgroup * r_time / gen_time;
            t->anc[t->root] = w_node2;
            t->time[w_node1] = r_time + outgroup * r_time;
            t->gen_length[w_node1] = t->time[w_node1] / gen_time;
            t->n_gen[w_node1] = t->gen_length[w_node1];
            t->n_replicas[w_node1] = 1;
            t->root = w_node2;
            t->n_leaves = leaves + 1;
            t->n_gleaves = leaves * ind_per_sp + 1;
        } else {
            t->n_leaves = leaves;
            t->n_gleaves = leaves * ind_per_sp;
        }
        t->locus_tree = NULL;
        free(i_nodes); free(node_ptrs);
        if (verbosity > 4) printf("\n\t\t\tDone");
    } else {
        /* ---- SSA ---- */
        if (verbosity > 3) printf("\n\tUsing SSA algorithm (fixed generations) to simulate the species tree...");
        while (iter <= MAX_IT) {
            ++iter;
            node_ptrs = sx_calloc((size_t)MAX_LEAVES, sizeof(int));
            t = sx_stree_new(1, 1, 0, 2, gen_time, Ne, mu); /* root = dummy node 0 */
            t->root = 0;
            *out_tree = t;
            w_node1 = sx_stree_add_node(t);
            SX_CHILD(t, t->root, 0) = w_node1;
            t->anc[w_node1] = t->root;
            anc_node = w_node1;
            extra_nodes = 0; c_time = 0;
            if (mrca_time == 0) {
                node_ptrs[0] = w_node1; n_leaves = 1; avail_leaves = 1; n_nodes = 1;
                t->n_gen[t->root] = 0;
                s_time = -log(gsl_rng_uniform_pos(seed)) / (w_prob * avail_leaves);
            } else {
                w_node1 = sx_stree_add_node(t); w_node2 = sx_stree_add_node(t);
                SX_CHILD(t, anc_node, 0) = w_node1; SX_CHILD(t, anc_node, 1) = w_node2;
                t->anc[w_node1] = anc_node; t->anc[w_node2] = anc_node;
                t->n_gen[anc_node] = 0; t->n_child[anc_node] = 2;
                n_leaves = 2; avail_leaves = 2; n_nodes = 3;
                node_ptrs[0] = w_node1; node_ptrs[1] = w_node2;
                s_time = -log(gsl_rng_uniform_pos(seed)) / (w_prob * avail_leaves);
            }
            while (c_time + s_time <= time && avail_leaves != 0) {
                if (n_leaves >= MAX_LEAVES) {
                    fprintf(stderr, "\nMaximum number of leaves reached. After checking your simulation parameters carefully (mainly species tree simulation birth rate) you can sort this problem out setting the environmental variable SIMPHY_MAXLEAVES with a value higher than %d\n", MAX_LEAVES);
                    return SX_MEM_ERROR;
                }
                if (avail_leaves > 1) node_index = (int)gsl_rng_uniform_int(seed, (unsigned long)avail_leaves);
                else node_index = 0;
                if (b_prob > gsl_rng_uniform(seed)) {
                    w_node1 = sx_stree_add_node(t); w_node2 = sx_stree_add_node(t);
                    anc_node = node_ptrs[node_index];
                    SX_CHILD(t, anc_node, 0) = w_node1; SX_CHILD(t, anc_node, 1) = w_node2;
                    t->anc[w_node1] = anc_node; t->anc[w_node2] = anc_node;
                    t->time[anc_node] = c_time + s_time;
                    t->n_gen[anc_node] = t->time[anc_node] / gen_time;
                    t->gen_length[anc_node] = t->n_gen[anc_node] - t->n_gen[t->anc[anc_node]];
                    t->n_child[anc_node] = 2;
                    node_ptrs[node_index] = w_node1; node_ptrs[avail_leaves] = w_node2;
                    ++avail_leaves; ++n_leaves; n_nodes += 2;
                } else {
                    w_node1 = node_ptrs[node_index];
                    if (avail_leaves == 1) { avail_leaves = 0; --n_leaves; n_nodes = 0; }
                    else {
                        anc_node = t->anc[w_node1];
                        if (SX_CHILD(t, anc_node, 0) == w_node1) w_node2 = SX_CHILD(t, anc_node, 1);
                        else w_node2 = SX_CHILD(t, anc_node, 0);
                        if (complete == 0) {
                            if (t->anc[anc_node] == t->root) avail_leaves = 0;
                            else {
                                t->anc[w_node2] = t->anc[anc_node];
                                if (SX_CHILD(t, t->anc[w_node2], 0) == anc_node) SX_CHILD(t, t->anc[w_node2], 0) = w_node2;
                                else SX_CHILD(t, t->anc[w_node2], 1) = w_node2;
                                if (t->anc[w_node2] >= 0) t->gen_length[w_node2] = t->n_gen[w_node2] - t->n_gen[t->anc[w_node2]];
                                node_ptrs[node_index] = node_ptrs[avail_leaves - 1];
                                /* w_node1 and anc_node become unreachable */
                                --avail_leaves; --n_leaves; n_nodes -= 2;
                            }
                        } else {
                            t->time[w_node1] = c_time + s_time;
                            t->n_gen[w_node1] = t->time[w_node1] / gen_time;
                            t->gen_length[w_node1] = t->n_gen[w_node1] - t->n_gen[anc_node];
                            node_ptrs[node_index] = node_ptrs[avail_leaves - 1];
                            --avail_leaves;
                            t->sp_index[w_node1] = n_leaves - avail_leaves;
                            n_nodes -= 2; extra_nodes += 2;
                        }
                    }
                }
                c_time += s_time;
                s_time = -log(gsl_rng_uniform_pos(seed)) / (w_prob * avail_leaves);
            }
            if (avail_leaves < 2) {
                sx_stree_free(t); t = NULL; *out_tree = NULL; free(node_ptrs); node_ptrs = NULL;
                if (verbosity > 2) printf("\n\t\t\tSpecies tree with less than 2 leaves or mrca smaller than expected, restart of the simulation of this tree. Try %d of %d", iter, MAX_IT);
                continue;
            } else { iter = 0; break; }
        }
        if (t == NULL) return SX_LOOP_ERROR;
        for (j = 0; j < avail_leaves; ++j) {
            w_node1 = node_ptrs[j];
            t->time[w_node1] = time;
            t->n_gen[w_node1] = time / gen_time;
            t->gen_length[w_node1] = t->n_gen[w_node1] - t->n_gen[t->anc[w_node1]];
            t->n_replicas[w_node1] = ind_per_sp;
            t->sp_index[w_node1] = j + n_leaves + 1 - avail_leaves;
        }
        if (outgroup > 0) {
            outgroup = time * outgroup;
            w_node1 = sx_stree_add_node(t);
            SX_CHILD(t, t->root, 1) = w_node1;
            t->anc[w_node1] = t->root;
            t->n_gen[t->root] = 0; t->n_child[t->root] = 2;
            t->gen_length[SX_CHILD(t, t->root, 0)] = outgroup / gen_time;
            t->time[w_node1] = time + outgroup;
            t->n_gen[w_node1] = t->time[w_node1] / gen_time;
            t->gen_length[w_node1] = t->n_gen[w_node1];
            t->n_replicas[w_node1] = 1;
            add_constant_ngen(t, SX_CHILD(t, t->root, 0), outgroup / gen_time, gen_time);
            t->n_leaves = avail_leaves + (extra_nodes / 2) + 1;
            t->n_gleaves = (avail_leaves + (extra_nodes / 2)) * ind_per_sp + 1;
        } else {
            w_node1 = SX_CHILD(t, t->root, 0);
            t->anc[w_node1] = -1;
            t->root = w_node1;   /* dummy root becomes unreachable */
            t->n_leaves = avail_leaves + (extra_nodes / 2);
            t->n_gleaves = (avail_leaves + (extra_nodes / 2) * ind_per_sp);
        }
        t->locus_tree = NULL;
        free(node_ptrs);
        if (iter > MAX_IT) return SX_LOOP_ERROR;
    }
    /* n_nodes is finalised by the collapse (unreachable arena nodes dropped). */
    return SX_NO_ERROR;
}

void sx_stree_rateheter(sx_stree *t, double alpha, gsl_rng *seed)
{
    for (int i = 0; i < t->n_nodes; ++i) t->mu_mult[i] *= gsl_ran_gamma(seed, alpha, 1 / alpha);
}
