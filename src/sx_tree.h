/* SPDX-License-Identifier: GPL-3.0-only */
/* Array-based (struct-of-arrays) species, locus and gene trees.
 * Node ids are positions in the arrays; children are stored with a fixed
 * stride of max_children (as in SimPhy) and -1 marks "no node". */
#ifndef SX_TREE_H
#define SX_TREE_H
#include "sx_common.h"

struct sx_ltree;

typedef struct sx_stree {
    int n_nodes, n_leaves, n_gleaves, max_children;
    double gen_time; int Ne; double mu;
    int root;
    int cap;
    int *index, *sp_index, *n_child, *n_replicas, *Ne_b, *anc, *child;
    double *n_gen, *time, *gen_length, *mu_mult, *gtime_mult;
    int *l_head, *n_lnodes;          /* locus-tree simulation linkage */
    struct sx_ltree *locus_tree;
    int n_nodes_hint;                /* verbosity>3: tree is printed by the caller after collapse */
} sx_stree;

typedef struct sx_ltree {
    int n_nodes, n_leaves, n_gleaves, max_children;
    double gen_time; int Ne; double mu;
    int root;
    int cap;
    int *index, *sp_index, *kind, *paralog, *n_child, *n_lin, *n_ilin, *n_olin, *fmax_nlin, *Ne_b, *anc, *lat, *conts, *child;
    double *n_gen, *time, *gen_length, *mu_mult, *gtime_mult;
    double **i_probs, **i_combprobs, **o_probs;
    sx_stree *species_tree;
} sx_ltree;

typedef struct sx_gtree {
    int n_nodes, max_children;
    double gen_time;
    int root;
    int *index, *sp_index, *replica, *paralog, *n_child, *anc, *contl, *conts, *child;
    double *n_gen, *gen_length, *bl;
    sx_ltree *locus_tree;
    sx_stree *species_tree;
} sx_gtree;

#define SX_CHILD(t, n, k) ((t)->child[(size_t)(n) * (size_t)(t)->max_children + (k)])

/* species tree */
sx_stree *sx_stree_new(int n_nodes, int n_leaves, int n_gleaves, int max_children, double gen_time, int Ne, double mu);
int sx_stree_add_node(sx_stree *t);
void sx_stree_free(sx_stree *t);
void sx_stree_collapse(sx_stree *t, int post_order);   /* CollapseSTree: reorder + index=i */
void sx_stree_reindex(sx_stree *t, int post_order);    /* ReindexSTree: labels only */
void sx_stree_reset_siml(sx_stree *t);                 /* ResetSTreeSimL */
long double sx_stree_gl_height(const sx_stree *t, int node);
long double sx_stree_cu_height(const sx_stree *t, int node, int g_Ne);
long double sx_stree_cu_length(const sx_stree *t, int node, int g_Ne);

/* locus tree */
sx_ltree *sx_ltree_new(int n_nodes, int n_leaves, int n_gleaves, int max_children, double gen_time, int Ne, double mu);
int sx_ltree_add_node(sx_ltree *t);
void sx_ltree_free(sx_ltree *t);
void sx_ltree_free_probs(sx_ltree *t);
/* Collapse reachable nodes from root in post-order into the first positions,
 * then set n_nodes = new_nnodes (>= reachable) and n_leaves. Unreachable
 * nodes are dropped; lat is cleared (relink==0 semantics). */
void sx_ltree_collapse(sx_ltree *t, int new_nnodes, int new_nleaves, int post_order);
void sx_ltree_copy_from_stree(sx_ltree *lt, sx_stree *st);   /* CopyStoLTree */

/* gene tree */
sx_gtree *sx_gtree_new(int n_nodes, int max_children, double gen_time);
void sx_gtree_free(sx_gtree *t);
void sx_gtree_reset(sx_gtree *t);                       /* MatchTrees reset part */
void sx_gtree_collapse(sx_gtree *t, int post_order);    /* CollapseGTree (relink handled by caller via pool) */
long double sx_gtree_cu_height(const sx_gtree *t, int node, int g_Ne);
long double sx_gtree_bl_height(const sx_gtree *t, int node);
long double sx_gtree_cu_length(const sx_gtree *t, int node, int g_Ne);
long double sx_gtree_bl_length(const sx_gtree *t, int node);

/* traversal order of the subtree rooted at root: order[k] = k-th visited node */
int sx_order_traversal(const int *n_child, const int *child, int max_children, int root, int post_order, int *order, int cap);

#endif

/* Deep copy of a species tree (arrays sized to n_nodes). */
sx_stree *sx_stree_clone(const sx_stree *src);
