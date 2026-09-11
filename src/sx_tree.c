#include "sx_tree.h"

/* ---------- traversal (iterative, children visited in stored order) ---------- */
int sx_order_traversal(const int *n_child, const int *child, int max_children, int root, int post_order, int *order, int cap)
{
    if (root < 0) return 0;
    int n = 0;
    /* stack of (node, next child index) */
    int scap = 64, sp = 0;
    int *st_node = sx_malloc(sizeof(int) * (size_t)scap), *st_k = sx_malloc(sizeof(int) * (size_t)scap);
    st_node[0] = root; st_k[0] = 0; sp = 1;
    if (!post_order) { if (n >= cap) sx_fatal("traversal overflow"); order[n++] = root; }
    while (sp > 0) {
        int v = st_node[sp - 1], k = st_k[sp - 1];
        if (k < n_child[v]) {
            st_k[sp - 1] = k + 1;
            int c = child[(size_t)v * (size_t)max_children + k];
            if (c < 0) sx_fatal("traversal: missing child");
            if (sp >= scap) { scap *= 2; st_node = sx_realloc(st_node, sizeof(int) * (size_t)scap); st_k = sx_realloc(st_k, sizeof(int) * (size_t)scap); }
            st_node[sp] = c; st_k[sp] = 0; ++sp;
            if (!post_order) { if (n >= cap) sx_fatal("traversal overflow"); order[n++] = c; }
        } else {
            --sp;
            if (post_order) { if (n >= cap) sx_fatal("traversal overflow"); order[n++] = v; }
        }
    }
    free(st_node); free(st_k);
    return n;
}

static void permute_int(int **arr, const int *order, int n, int cap)
{
    int *nw = sx_calloc((size_t)cap, sizeof(int));
    for (int k = 0; k < n; ++k) nw[k] = (*arr)[order[k]];
    free(*arr); *arr = nw;
}
static void permute_dbl(double **arr, const int *order, int n, int cap)
{
    double *nw = sx_calloc((size_t)cap, sizeof(double));
    for (int k = 0; k < n; ++k) nw[k] = (*arr)[order[k]];
    free(*arr); *arr = nw;
}
static void permute_ptr(double ***arr, const int *order, int n, int cap)
{
    double **nw = sx_calloc((size_t)cap, sizeof(double *));
    for (int k = 0; k < n; ++k) nw[k] = (*arr)[order[k]];
    free(*arr); *arr = nw;
}
static void permute_child(int **child, const int *order, const int *newpos, const int *n_child_old, int n, int cap, int mc)
{
    int *nw = sx_malloc(sizeof(int) * (size_t)cap * (size_t)mc);
    for (size_t i = 0; i < (size_t)cap * (size_t)mc; ++i) nw[i] = -1;
    for (int k = 0; k < n; ++k) {
        int old = order[k];
        for (int j = 0; j < n_child_old[old]; ++j) {
            int c = (*child)[(size_t)old * (size_t)mc + j];
            nw[(size_t)k * (size_t)mc + j] = c >= 0 ? newpos[c] : -1;
        }
    }
    free(*child); *child = nw;
}
static void remap_ids(int *arr, const int *newpos, int n)
{
    for (int k = 0; k < n; ++k) arr[k] = arr[k] >= 0 ? newpos[arr[k]] : -1;
}

/* ---------- species tree ---------- */
static void stree_init_node(sx_stree *t, int i)
{
    t->index[i] = i; t->sp_index[i] = 0; t->n_child[i] = 0; t->n_replicas[i] = 0; t->Ne_b[i] = 0; t->anc[i] = -1;
    for (int j = 0; j < t->max_children; ++j) SX_CHILD(t, i, j) = -1;
    t->n_gen[i] = 0; t->time[i] = 0; t->gen_length[i] = 0; t->mu_mult[i] = 1; t->gtime_mult[i] = 1;
    t->l_head[i] = -1; t->n_lnodes[i] = 0;
}

static void stree_alloc(sx_stree *t, int cap)
{
    t->cap = cap;
    t->index = sx_calloc((size_t)cap, sizeof(int)); t->sp_index = sx_calloc((size_t)cap, sizeof(int));
    t->n_child = sx_calloc((size_t)cap, sizeof(int)); t->n_replicas = sx_calloc((size_t)cap, sizeof(int));
    t->Ne_b = sx_calloc((size_t)cap, sizeof(int)); t->anc = sx_calloc((size_t)cap, sizeof(int));
    t->child = sx_calloc((size_t)cap * (size_t)t->max_children, sizeof(int));
    t->n_gen = sx_calloc((size_t)cap, sizeof(double)); t->time = sx_calloc((size_t)cap, sizeof(double));
    t->gen_length = sx_calloc((size_t)cap, sizeof(double)); t->mu_mult = sx_calloc((size_t)cap, sizeof(double));
    t->gtime_mult = sx_calloc((size_t)cap, sizeof(double));
    t->l_head = sx_calloc((size_t)cap, sizeof(int)); t->n_lnodes = sx_calloc((size_t)cap, sizeof(int));
}

sx_stree *sx_stree_new(int n_nodes, int n_leaves, int n_gleaves, int max_children, double gen_time, int Ne, double mu)
{
    sx_stree *t = sx_calloc(1, sizeof(*t));
    t->n_nodes = n_nodes; t->n_leaves = n_leaves; t->n_gleaves = n_gleaves; t->max_children = max_children;
    t->gen_time = gen_time; t->Ne = Ne; t->mu = mu; t->root = -1; t->locus_tree = NULL;
    stree_alloc(t, n_nodes > 0 ? n_nodes : 16);
    for (int i = 0; i < n_nodes; ++i) stree_init_node(t, i);
    return t;
}

static void stree_grow(sx_stree *t, int need)
{
    if (need <= t->cap) return;
    int nc = t->cap * 2; if (nc < need) nc = need;
    t->index = sx_realloc(t->index, sizeof(int) * (size_t)nc); t->sp_index = sx_realloc(t->sp_index, sizeof(int) * (size_t)nc);
    t->n_child = sx_realloc(t->n_child, sizeof(int) * (size_t)nc); t->n_replicas = sx_realloc(t->n_replicas, sizeof(int) * (size_t)nc);
    t->Ne_b = sx_realloc(t->Ne_b, sizeof(int) * (size_t)nc); t->anc = sx_realloc(t->anc, sizeof(int) * (size_t)nc);
    t->child = sx_realloc(t->child, sizeof(int) * (size_t)nc * (size_t)t->max_children);
    t->n_gen = sx_realloc(t->n_gen, sizeof(double) * (size_t)nc); t->time = sx_realloc(t->time, sizeof(double) * (size_t)nc);
    t->gen_length = sx_realloc(t->gen_length, sizeof(double) * (size_t)nc); t->mu_mult = sx_realloc(t->mu_mult, sizeof(double) * (size_t)nc);
    t->gtime_mult = sx_realloc(t->gtime_mult, sizeof(double) * (size_t)nc);
    t->l_head = sx_realloc(t->l_head, sizeof(int) * (size_t)nc); t->n_lnodes = sx_realloc(t->n_lnodes, sizeof(int) * (size_t)nc);
    t->cap = nc;
}

int sx_stree_add_node(sx_stree *t)
{
    stree_grow(t, t->n_nodes + 1);
    int i = t->n_nodes++;
    stree_init_node(t, i);
    return i;
}

void sx_stree_free(sx_stree *t)
{
    if (!t) return;
    free(t->index); free(t->sp_index); free(t->n_child); free(t->n_replicas); free(t->Ne_b); free(t->anc); free(t->child);
    free(t->n_gen); free(t->time); free(t->gen_length); free(t->mu_mult); free(t->gtime_mult); free(t->l_head); free(t->n_lnodes);
    free(t);
}

void sx_stree_collapse(sx_stree *t, int post_order)
{
    if (t->root < 0) sx_fatal("CollapseSTree: no root");
    int *order = sx_malloc(sizeof(int) * (size_t)t->n_nodes);
    int n = sx_order_traversal(t->n_child, t->child, t->max_children, t->root, post_order, order, t->n_nodes);
    int *newpos = sx_malloc(sizeof(int) * (size_t)t->n_nodes);
    for (int i = 0; i < t->n_nodes; ++i) newpos[i] = -1;
    for (int k = 0; k < n; ++k) newpos[order[k]] = k;
    int cap = t->cap;
    permute_child(&t->child, order, newpos, t->n_child, n, cap, t->max_children);
    permute_int(&t->sp_index, order, n, cap); permute_int(&t->n_child, order, n, cap); permute_int(&t->n_replicas, order, n, cap);
    permute_int(&t->Ne_b, order, n, cap); permute_int(&t->anc, order, n, cap); remap_ids(t->anc, newpos, n);
    permute_dbl(&t->n_gen, order, n, cap); permute_dbl(&t->time, order, n, cap); permute_dbl(&t->gen_length, order, n, cap);
    permute_dbl(&t->mu_mult, order, n, cap); permute_dbl(&t->gtime_mult, order, n, cap);
    permute_int(&t->l_head, order, n, cap); permute_int(&t->n_lnodes, order, n, cap);
    for (int k = 0; k < n; ++k) t->index[k] = k;
    t->root = newpos[t->root];
    t->n_nodes = n;
    free(order); free(newpos);
}

void sx_stree_reindex(sx_stree *t, int post_order)
{
    int *order = sx_malloc(sizeof(int) * (size_t)t->n_nodes);
    int n = sx_order_traversal(t->n_child, t->child, t->max_children, t->root, post_order, order, t->n_nodes);
    for (int k = 0; k < n; ++k) t->index[order[k]] = k;
    free(order);
}

void sx_stree_reset_siml(sx_stree *t)
{
    for (int i = 0; i < t->n_nodes; ++i) { t->index[i] = i; t->l_head[i] = -1; t->n_lnodes[i] = 0; }
    t->locus_tree = NULL;
}

long double sx_stree_gl_height(const sx_stree *t, int node)
{
    long double current = 0, nw;
    for (int i = 0; i < t->n_child[node]; ++i) { nw = sx_stree_gl_height(t, SX_CHILD(t, node, i)); if (nw > current) current = nw; }
    return current + t->gen_length[node];
}
long double sx_stree_cu_height(const sx_stree *t, int node, int g_Ne)
{
    long double current = 0, nw;
    for (int i = 0; i < t->n_child[node]; ++i) { nw = sx_stree_cu_height(t, SX_CHILD(t, node, i), g_Ne); if (nw > current) current = nw; }
    return t->gen_length[node] / (t->Ne_b[node] != 0 ? t->Ne_b[node] : g_Ne) + current;
}
long double sx_stree_cu_length(const sx_stree *t, int node, int g_Ne)
{
    long double sum = 0;
    for (int i = 0; i < t->n_child[node]; ++i) sum += sx_stree_cu_length(t, SX_CHILD(t, node, i), g_Ne);
    return t->gen_length[node] / (t->Ne_b[node] != 0 ? t->Ne_b[node] : g_Ne) + sum;
}

/* ---------- locus tree ---------- */
static void ltree_init_node(sx_ltree *t, int i)
{
    t->index[i] = i; t->sp_index[i] = 0; t->kind[i] = SP; t->paralog[i] = 0; t->n_child[i] = 0; t->n_lin[i] = 0;
    t->n_ilin[i] = 0; t->n_olin[i] = 0; t->fmax_nlin[i] = 0; t->Ne_b[i] = 0; t->anc[i] = -1; t->lat[i] = -1; t->conts[i] = -1;
    for (int j = 0; j < t->max_children; ++j) SX_CHILD(t, i, j) = -1;
    t->n_gen[i] = 0; t->time[i] = 0; t->gen_length[i] = 0; t->mu_mult[i] = 1; t->gtime_mult[i] = 1;
    t->i_probs[i] = NULL; t->i_combprobs[i] = NULL; t->o_probs[i] = NULL;
}
static void ltree_alloc(sx_ltree *t, int cap)
{
    t->cap = cap;
#define A(f) t->f = sx_calloc((size_t)cap, sizeof(*t->f))
    A(index); A(sp_index); A(kind); A(paralog); A(n_child); A(n_lin); A(n_ilin); A(n_olin); A(fmax_nlin); A(Ne_b); A(anc); A(lat); A(conts);
    A(n_gen); A(time); A(gen_length); A(mu_mult); A(gtime_mult); A(i_probs); A(i_combprobs); A(o_probs);
#undef A
    t->child = sx_calloc((size_t)cap * (size_t)t->max_children, sizeof(int));
}
sx_ltree *sx_ltree_new(int n_nodes, int n_leaves, int n_gleaves, int max_children, double gen_time, int Ne, double mu)
{
    sx_ltree *t = sx_calloc(1, sizeof(*t));
    t->n_nodes = n_nodes; t->n_leaves = n_leaves; t->n_gleaves = n_gleaves; t->max_children = max_children;
    t->gen_time = gen_time; t->Ne = Ne; t->mu = mu; t->root = -1; t->species_tree = NULL;
    ltree_alloc(t, n_nodes > 0 ? n_nodes : 16);
    for (int i = 0; i < n_nodes; ++i) ltree_init_node(t, i);
    if (n_nodes == 1) t->root = 0;
    return t;
}
static void ltree_grow(sx_ltree *t, int need)
{
    if (need <= t->cap) return;
    int nc = t->cap * 2; if (nc < need) nc = need;
#define R(f) t->f = sx_realloc(t->f, sizeof(*t->f) * (size_t)nc)
    R(index); R(sp_index); R(kind); R(paralog); R(n_child); R(n_lin); R(n_ilin); R(n_olin); R(fmax_nlin); R(Ne_b); R(anc); R(lat); R(conts);
    R(n_gen); R(time); R(gen_length); R(mu_mult); R(gtime_mult); R(i_probs); R(i_combprobs); R(o_probs);
#undef R
    t->child = sx_realloc(t->child, sizeof(int) * (size_t)nc * (size_t)t->max_children);
    t->cap = nc;
}
int sx_ltree_add_node(sx_ltree *t)
{
    ltree_grow(t, t->n_nodes + 1);
    int i = t->n_nodes++;
    ltree_init_node(t, i);
    return i;
}
void sx_ltree_free_probs(sx_ltree *t)
{
    for (int i = 0; i < t->n_nodes; ++i) {
        free(t->i_probs[i]); free(t->i_combprobs[i]); free(t->o_probs[i]);
        t->i_probs[i] = t->i_combprobs[i] = t->o_probs[i] = NULL;
    }
}
void sx_ltree_free(sx_ltree *t)
{
    if (!t) return;
    sx_ltree_free_probs(t);
    free(t->index); free(t->sp_index); free(t->kind); free(t->paralog); free(t->n_child); free(t->n_lin); free(t->n_ilin); free(t->n_olin);
    free(t->fmax_nlin); free(t->Ne_b); free(t->anc); free(t->lat); free(t->conts); free(t->child);
    free(t->n_gen); free(t->time); free(t->gen_length); free(t->mu_mult); free(t->gtime_mult);
    free(t->i_probs); free(t->i_combprobs); free(t->o_probs);
    free(t);
}

void sx_ltree_collapse(sx_ltree *t, int new_nnodes, int new_nleaves, int post_order)
{
    if (t->root < 0) sx_fatal("CollapseLTree: no root");
    int *order = sx_malloc(sizeof(int) * (size_t)t->n_nodes);
    int n = sx_order_traversal(t->n_child, t->child, t->max_children, t->root, post_order, order, t->n_nodes);
    if (new_nnodes < n) sx_fatal("CollapseResizeLTree: new size %d smaller than reachable nodes %d", new_nnodes, n);
    int *newpos = sx_malloc(sizeof(int) * (size_t)t->n_nodes);
    for (int i = 0; i < t->n_nodes; ++i) newpos[i] = -1;
    for (int k = 0; k < n; ++k) newpos[order[k]] = k;
    /* free probs of unreachable nodes */
    for (int i = 0; i < t->n_nodes; ++i) if (newpos[i] < 0) { free(t->i_probs[i]); free(t->i_combprobs[i]); free(t->o_probs[i]); t->i_probs[i] = t->i_combprobs[i] = t->o_probs[i] = NULL; }
    int cap = t->cap > new_nnodes ? t->cap : new_nnodes;
    permute_child(&t->child, order, newpos, t->n_child, n, cap, t->max_children);
#define P(f) permute_int(&t->f, order, n, cap)
    P(sp_index); P(kind); P(paralog); P(n_child); P(n_lin); P(n_ilin); P(n_olin); P(fmax_nlin); P(Ne_b); P(anc); P(lat); P(conts);
#undef P
    remap_ids(t->anc, newpos, n);
    for (int k = 0; k < n; ++k) t->lat[k] = -1;   /* relink == 0 */
    permute_dbl(&t->n_gen, order, n, cap); permute_dbl(&t->time, order, n, cap); permute_dbl(&t->gen_length, order, n, cap);
    permute_dbl(&t->mu_mult, order, n, cap); permute_dbl(&t->gtime_mult, order, n, cap);
    permute_ptr(&t->i_probs, order, n, cap); permute_ptr(&t->i_combprobs, order, n, cap); permute_ptr(&t->o_probs, order, n, cap);
    t->index = sx_realloc(t->index, sizeof(int) * (size_t)cap);
    t->cap = cap;
    for (int k = 0; k < n; ++k) t->index[k] = k;
    t->root = newpos[t->root];
    t->n_nodes = n;
    for (int k = n; k < new_nnodes; ++k) { t->n_nodes = k + 1; ltree_init_node(t, k); }
    t->n_nodes = new_nnodes;
    t->n_leaves = new_nleaves;
    free(order); free(newpos);
}

/* CopyStoLTree: direct conversion of a (pre-ordered) species tree into a
 * post-ordered locus tree.  The species tree indices are relabelled to
 * post-order during the copy and restored to pre-order afterwards, as in SimPhy. */
void sx_ltree_copy_from_stree(sx_ltree *lt, sx_stree *st)
{
    if (st->max_children != lt->max_children || st->n_gleaves != lt->n_gleaves || st->n_leaves != lt->n_leaves || st->n_nodes != lt->n_nodes)
        sx_check(SX_MEM_ERROR, "CopyStoLTree: tree size mismatch");
    int *order = sx_malloc(sizeof(int) * (size_t)st->n_nodes);
    int n = sx_order_traversal(st->n_child, st->child, st->max_children, st->root, 1, order, st->n_nodes);
    int *post = sx_malloc(sizeof(int) * (size_t)st->n_nodes);
    for (int k = 0; k < n; ++k) post[order[k]] = k;
    for (int k = 0; k < n; ++k) {
        int s = order[k];
        lt->sp_index[k] = st->sp_index[s]; lt->n_child[k] = st->n_child[s]; lt->n_lin[k] = st->n_replicas[s];
        lt->Ne_b[k] = st->Ne_b[s]; lt->n_gen[k] = st->n_gen[s]; lt->time[k] = st->time[s]; lt->gen_length[k] = st->gen_length[s];
        lt->conts[k] = s; lt->mu_mult[k] = st->mu_mult[s]; lt->gtime_mult[k] = st->gtime_mult[s];
        st->n_lnodes[s] = 1; st->l_head[s] = k;
        for (int j = 0; j < st->n_child[s]; ++j) SX_CHILD(lt, k, j) = post[SX_CHILD(st, s, j)];
        lt->anc[k] = st->anc[s] >= 0 ? post[st->anc[s]] : -1;
    }
    lt->root = st->l_head[st->root];
    lt->gen_time = st->gen_time; lt->Ne = st->Ne; lt->mu = st->mu;
    lt->species_tree = st; st->locus_tree = lt;
    free(order); free(post);
}

/* ---------- gene tree ---------- */
sx_gtree *sx_gtree_new(int n_nodes, int max_children, double gen_time)
{
    sx_gtree *t = sx_calloc(1, sizeof(*t));
    t->n_nodes = n_nodes; t->max_children = max_children; t->gen_time = gen_time; t->root = -1;
    size_t n = (size_t)(n_nodes > 0 ? n_nodes : 1);
    t->index = sx_calloc(n, sizeof(int)); t->sp_index = sx_calloc(n, sizeof(int)); t->replica = sx_calloc(n, sizeof(int));
    t->paralog = sx_calloc(n, sizeof(int)); t->n_child = sx_calloc(n, sizeof(int)); t->anc = sx_calloc(n, sizeof(int));
    t->contl = sx_calloc(n, sizeof(int)); t->conts = sx_calloc(n, sizeof(int));
    t->child = sx_calloc(n * (size_t)max_children, sizeof(int));
    t->n_gen = sx_calloc(n, sizeof(double)); t->gen_length = sx_calloc(n, sizeof(double)); t->bl = sx_calloc(n, sizeof(double));
    sx_gtree_reset(t);
    return t;
}
void sx_gtree_reset(sx_gtree *t)
{
    for (int i = 0; i < t->n_nodes; ++i) {
        t->index[i] = i; t->sp_index[i] = 0; t->paralog[i] = 0; t->replica[i] = 0; t->conts[i] = -1; t->contl[i] = -1;
        t->n_child[i] = 0; t->anc[i] = -1; t->n_gen[i] = 0.0; t->bl[i] = 0.0; t->gen_length[i] = 0.0;
        for (int j = 0; j < t->max_children; ++j) SX_CHILD(t, i, j) = -1;
    }
    t->root = -1;
}
void sx_gtree_free(sx_gtree *t)
{
    if (!t) return;
    free(t->index); free(t->sp_index); free(t->replica); free(t->paralog); free(t->n_child); free(t->anc); free(t->contl); free(t->conts);
    free(t->child); free(t->n_gen); free(t->gen_length); free(t->bl);
    free(t);
}
void sx_gtree_collapse(sx_gtree *t, int post_order)
{
    if (t->root < 0) sx_fatal("CollapseGTree: no root");
    int *order = sx_malloc(sizeof(int) * (size_t)t->n_nodes);
    int n = sx_order_traversal(t->n_child, t->child, t->max_children, t->root, post_order, order, t->n_nodes);
    if (n != t->n_nodes) sx_fatal("CollapseGTree: reachable nodes %d != n_nodes %d", n, t->n_nodes);
    int *newpos = sx_malloc(sizeof(int) * (size_t)n);
    for (int k = 0; k < n; ++k) newpos[order[k]] = k;
    permute_child(&t->child, order, newpos, t->n_child, n, n, t->max_children);
    permute_int(&t->sp_index, order, n, n); permute_int(&t->replica, order, n, n); permute_int(&t->paralog, order, n, n);
    permute_int(&t->n_child, order, n, n); permute_int(&t->anc, order, n, n); remap_ids(t->anc, newpos, n);
    permute_int(&t->contl, order, n, n); permute_int(&t->conts, order, n, n);
    permute_dbl(&t->n_gen, order, n, n); permute_dbl(&t->gen_length, order, n, n); permute_dbl(&t->bl, order, n, n);
    for (int k = 0; k < n; ++k) t->index[k] = k;
    t->root = newpos[t->root];
    free(order); free(newpos);
}
long double sx_gtree_cu_height(const sx_gtree *t, int node, int g_Ne)
{
    long double current = 0, nw;
    for (int i = 0; i < t->n_child[node]; ++i) { nw = sx_gtree_cu_height(t, SX_CHILD(t, node, i), g_Ne); if (nw > current) current = nw; }
    int l = t->contl[node];
    int Ne = t->locus_tree->Ne_b[l];
    return t->gen_length[node] / (Ne != 0 ? Ne : g_Ne) + current;
}
long double sx_gtree_bl_height(const sx_gtree *t, int node)
{
    long double current = 0, nw;
    for (int i = 0; i < t->n_child[node]; ++i) { nw = sx_gtree_bl_height(t, SX_CHILD(t, node, i)); if (nw > current) current = nw; }
    return current + t->bl[node];
}
long double sx_gtree_cu_length(const sx_gtree *t, int node, int g_Ne)
{
    long double sum = 0;
    for (int i = 0; i < t->n_child[node]; ++i) sum += sx_gtree_cu_length(t, SX_CHILD(t, node, i), g_Ne);
    int l = t->contl[node];
    int Ne = t->locus_tree->Ne_b[l];
    return t->gen_length[node] / (Ne != 0 ? Ne : g_Ne) + sum;
}
long double sx_gtree_bl_length(const sx_gtree *t, int node)
{
    long double sum = 0;
    for (int i = 0; i < t->n_child[node]; ++i) sum += sx_gtree_bl_length(t, SX_CHILD(t, node, i));
    return sum + t->bl[node];
}

sx_stree *sx_stree_clone(const sx_stree *s)
{
    sx_stree *t = sx_stree_new(s->n_nodes, s->n_leaves, s->n_gleaves, s->max_children, s->gen_time, s->Ne, s->mu);
    size_t n = (size_t)s->n_nodes;
    memcpy(t->index, s->index, n * sizeof(int)); memcpy(t->sp_index, s->sp_index, n * sizeof(int)); memcpy(t->n_child, s->n_child, n * sizeof(int));
    memcpy(t->n_replicas, s->n_replicas, n * sizeof(int)); memcpy(t->Ne_b, s->Ne_b, n * sizeof(int)); memcpy(t->anc, s->anc, n * sizeof(int));
    memcpy(t->child, s->child, n * (size_t)s->max_children * sizeof(int));
    memcpy(t->n_gen, s->n_gen, n * sizeof(double)); memcpy(t->time, s->time, n * sizeof(double)); memcpy(t->gen_length, s->gen_length, n * sizeof(double));
    memcpy(t->mu_mult, s->mu_mult, n * sizeof(double)); memcpy(t->gtime_mult, s->gtime_mult, n * sizeof(double));
    memcpy(t->l_head, s->l_head, n * sizeof(int)); memcpy(t->n_lnodes, s->n_lnodes, n * sizeof(int));
    t->root = s->root; t->locus_tree = NULL;
    return t;
}
