#include "sx_write.h"

static inline const char *nm(const sx_names *names, int idx) { return names->names[idx]; }

/* Generic iterative Newick emission. leaf(b,node) prints a leaf; inode(b,node)
 * prints the closing of an internal node (excluding the root's ");\n"). */
typedef struct {
    sx_buf *b;
    const void *tree;
    const sx_names *names;
    int time_units, int_labels;
} wctx;

static void emit_newick(sx_buf *b, int root, const int *n_child, const int *child, int mc, const int *anc,
                        void (*leaf)(wctx *, int), void (*inode_close)(wctx *, int), wctx *cx)
{
    int scap = 64, sp = 0;
    int *st_node = sx_malloc(sizeof(int) * (size_t)scap), *st_k = sx_malloc(sizeof(int) * (size_t)scap);
    st_node[0] = root; st_k[0] = 0; sp = 1;
    if (n_child[root] == 0) { leaf(cx, root); }
    else sx_buf_putc(b, '(');
    while (sp > 0) {
        int v = st_node[sp - 1], k = st_k[sp - 1];
        if (n_child[v] == 0) { --sp; continue; }
        if (k < n_child[v]) {
            if (k > 0) sx_buf_putc(b, ',');
            st_k[sp - 1] = k + 1;
            int c = child[(size_t)v * (size_t)mc + k];
            if (n_child[c] == 0) { leaf(cx, c); }
            else {
                sx_buf_putc(b, '(');
                if (sp >= scap) { scap *= 2; st_node = sx_realloc(st_node, sizeof(int) * (size_t)scap); st_k = sx_realloc(st_k, sizeof(int) * (size_t)scap); }
                st_node[sp] = c; st_k[sp] = 0; ++sp;
            }
        } else {
            --sp;
            if (anc[v] != -1) inode_close(cx, v);
            else sx_buf_puts(b, ");\n");
        }
    }
    free(st_node); free(st_k);
}

/* ---- species tree ---- */
static void s_leaf(wctx *cx, int p)
{
    const sx_stree *t = cx->tree;
    double len = cx->time_units ? t->gen_length[p] * t->gtime_mult[p] * t->gen_time : t->gen_length[p];
    if (cx->names == NULL) sx_buf_printf(cx->b, "%d:%.8lf", t->sp_index[p], len);
    else sx_buf_printf(cx->b, "%s:%.8lf", nm(cx->names, t->sp_index[p]), len);
}
static void s_close(wctx *cx, int p)
{
    const sx_stree *t = cx->tree;
    double len = cx->time_units ? t->gen_length[p] * t->gtime_mult[p] * t->gen_time : t->gen_length[p];
    if (cx->int_labels) sx_buf_printf(cx->b, ")%d:%.8lf", t->index[p], len);
    else sx_buf_printf(cx->b, "):%.8lf", len);
}
void sx_write_stree(sx_buf *b, const sx_stree *t, const sx_names *names, int time_units, int int_labels)
{
    if (t->root < 0) sx_check(SX_MEM_ERROR, "WriteSTree");
    if (time_units && t->gen_time == 0) sx_check(SX_UNEXPECTED_VALUE, "WriteSTree");
    wctx cx = { b, t, names, time_units, int_labels };
    emit_newick(b, t->root, t->n_child, t->child, t->max_children, t->anc, s_leaf, s_close, &cx);
}

/* ---- locus tree ---- */
static void l_leaf(wctx *cx, int p)
{
    const sx_ltree *t = cx->tree;
    const sx_stree *st = t->species_tree;
    double len = cx->time_units ? t->gen_length[p] * t->gen_time * t->gtime_mult[p] : t->gen_length[p];
    const char *pre = NULL;
    switch (t->kind[p]) {
    default:
        if (cx->names == NULL) sx_buf_printf(cx->b, "%d_%d:%.8lf", t->sp_index[p], t->paralog[p], len);
        else sx_buf_printf(cx->b, "%s_%d:%.8lf", nm(cx->names, t->sp_index[p]), t->paralog[p], len);
        return;
    case LOSS: pre = "Lost"; break;
    case RTRFR: pre = "Rtransf"; break;
    case RGC: pre = "Rgc"; break;
    }
    if (t->conts[p] == -1) sx_buf_printf(cx->b, "%s-%d_%d:%.8lf", pre, t->index[p], t->paralog[p], len);
    else if (st->n_child[t->conts[p]] == 0 && cx->names != NULL) sx_buf_printf(cx->b, "%s-%s_%d:%.8lf", pre, nm(cx->names, t->sp_index[p]), t->paralog[p], len);
    else sx_buf_printf(cx->b, "%s-%d_%d:%.8lf", pre, st->index[t->conts[p]], t->paralog[p], len);
}
static void l_close(wctx *cx, int p)
{
    const sx_ltree *t = cx->tree;
    double len = cx->time_units ? t->gen_length[p] * t->gen_time * t->gtime_mult[p] : t->gen_length[p];
    if (cx->int_labels) sx_buf_printf(cx->b, ")%d:%.8lf", t->index[p], len);
    else sx_buf_printf(cx->b, "):%.8lf", len);
}
void sx_write_ltree(sx_buf *b, const sx_ltree *t, const sx_names *names, int time_units, int int_labels)
{
    if (t->root < 0) sx_check(SX_MEM_ERROR, "WriteLTree");
    if (time_units && t->gen_time == 0) sx_check(SX_UNEXPECTED_VALUE, "WriteLTree");
    wctx cx = { b, t, names, time_units, int_labels };
    emit_newick(b, t->root, t->n_child, t->child, t->max_children, t->anc, l_leaf, l_close, &cx);
}

/* ---- gene tree ---- */
static void g_leaf(wctx *cx, int p)
{
    const sx_gtree *t = cx->tree;
    const sx_ltree *lt = t->locus_tree;
    const sx_stree *st = t->species_tree;
    const char *pre = NULL;
    switch (lt->kind[t->contl[p]]) {
    default:
        if (cx->names == NULL) sx_buf_printf(cx->b, "%d_%d_%d:%.8lf", t->sp_index[p], t->paralog[p], t->replica[p], t->bl[p]);
        else sx_buf_printf(cx->b, "%s_%d_%d:%.8lf", nm(cx->names, t->sp_index[p]), t->paralog[p], t->replica[p], t->bl[p]);
        return;
    case LOSS: pre = "Lost"; break;
    case RTRFR: pre = "Rtransf"; break;
    case RGC: pre = "Rgc"; break;
    }
    if (st->n_child[t->conts[p]] == 0 && cx->names != NULL) sx_buf_printf(cx->b, "%s-%s_%d_0:%.8lf", pre, nm(cx->names, t->sp_index[p]), t->paralog[p], t->bl[p]);
    else sx_buf_printf(cx->b, "%s-%d_%d_0:%.8lf", pre, st->index[t->conts[p]], t->paralog[p], t->bl[p]);
}
static void g_close(wctx *cx, int p)
{
    const sx_gtree *t = cx->tree;
    if (cx->int_labels) sx_buf_printf(cx->b, ")%d:%.8lf", t->index[p], t->bl[p]);
    else sx_buf_printf(cx->b, "):%.8lf", t->bl[p]);
}
void sx_write_gtree(sx_buf *b, const sx_gtree *t, const sx_names *names, int int_labels)
{
    if (t->root < 0) sx_check(SX_MEM_ERROR, "WriteGTree");
    wctx cx = { b, t, names, 0, int_labels };
    emit_newick(b, t->root, t->n_child, t->child, t->max_children, t->anc, g_leaf, g_close, &cx);
}

/* ---- daughters ---- */
void sx_write_daughters(sx_buf *b, const sx_ltree *t, const sx_names *names)
{
    if (t->root < 0) sx_check(SX_MEM_ERROR, "WriteDaughters");
    int *order = sx_malloc(sizeof(int) * (size_t)t->n_nodes);
    int n = sx_order_traversal(t->n_child, t->child, t->max_children, t->root, 1, order, t->n_nodes);
    size_t start = b->len;
    for (int k = 0; k < n; ++k) {
        int p = order[k];
        if (t->n_child[p] == 0) continue;
        switch (t->kind[p]) {
        case DUP: case TRFR: case GC: {
            int d = SX_CHILD(t, p, 1);
            if (t->n_child[d] != 0) sx_buf_printf(b, "%d,", t->index[d]);
            else if (names != NULL) sx_buf_printf(b, "'%s_%d',", nm(names, t->sp_index[d]), t->paralog[d]);
            else sx_buf_printf(b, "'%d_%d',", t->species_tree->index[t->conts[d]], t->paralog[d]);
            break;
        }
        default: break;
        }
    }
    free(order);
    /* fseek(-1, SEEK_END); fprintf("\n") */
    if (b->len > start) b->data[b->len - 1] = '\n';
    else if (b->len > 0) b->data[b->len - 1] = '\n';
    else sx_buf_putc(b, '\n');
}

/* ---- mappings ---- */
long sx_write_mapping_sl(sx_buf *b, sx_stree *st, const sx_ltree *lt, const sx_names *names)
{
    sx_buf_puts(b, "Lt_node\tLt_paralog\tLt_event\tSt_node\n");
    sx_stree_reindex(st, 1);
    for (int i = 0; i < lt->n_nodes; ++i) {
        int c = lt->conts[i];
        if (c == -1) return SX_UNEXPECTED_VALUE;
        if (lt->n_child[i] != 0) {
            const char *ev;
            switch (lt->kind[i]) {
            case SP: ev = "Sp"; break; case DUP: ev = "Dup"; break; case TRFR: ev = "Transf"; break; case GC: ev = "Gconv"; break;
            default: return SX_UNEXPECTED_VALUE;
            }
            if (st->n_child[c] != 0) sx_buf_printf(b, "%d\t%d\t%s\t%d\n", lt->index[i], lt->paralog[i], ev, st->index[c]);
            else if (names == NULL) sx_buf_printf(b, "%d\t%d\t%s\t'%d'\n", lt->index[i], lt->paralog[i], ev, lt->sp_index[i]);
            else sx_buf_printf(b, "%d\t%d\t%s\t'%s'\n", lt->index[i], lt->paralog[i], ev, nm(names, lt->sp_index[i]));
        } else if (lt->kind[i] == SP) {
            if (names == NULL) sx_buf_printf(b, "'%d_%d'\t%d\tLeaf\t'%d'\n", lt->sp_index[i], lt->paralog[i], lt->paralog[i], lt->sp_index[i]);
            else sx_buf_printf(b, "'%s_%d'\t%d\tLeaf\t'%s'\n", nm(names, lt->sp_index[i]), lt->paralog[i], lt->paralog[i], nm(names, lt->sp_index[i]));
        } else if (lt->kind[i] == LOSS) sx_buf_printf(b, "'Lost-%d_%d'\t%d\tLoss\t%d\n", st->index[c], lt->paralog[i], lt->paralog[i], st->index[c]);
        else if (lt->kind[i] == RTRFR) sx_buf_printf(b, "'RTransf-%d_%d'\t%d\tRTransf\t%d\n", st->index[c], lt->paralog[i], lt->paralog[i], st->index[c]);
        else sx_buf_printf(b, "'Rgc-%d_%d'\t%d\tRgc\t%d\n", st->index[c], lt->paralog[i], lt->paralog[i], st->index[c]);
    }
    return SX_NO_ERROR;
}

long sx_write_mapping_lg(sx_buf *b, const sx_gtree *gt, const sx_names *names)
{
    const sx_ltree *lt = gt->locus_tree;
    const sx_stree *st = gt->species_tree;
    sx_buf_puts(b, "Gt_node\tLt_node\tLt_paralog\n");
    for (int i = 0; i < gt->n_nodes; ++i) {
        int l = gt->contl[i];
        if (l == -1) return SX_UNEXPECTED_VALUE;
        if (gt->n_child[i] != 0) {
            if (lt->n_child[l] != 0) sx_buf_printf(b, "%d\t%d\t%d\n", gt->index[i], lt->index[l], gt->paralog[i]);
            else {
                switch (lt->kind[l]) {
                case LOSS: sx_buf_printf(b, "%d\t'Lost-%d_%d'\t%d\n", gt->index[i], st->index[gt->conts[i]], gt->paralog[i], gt->paralog[i]); break;
                case RTRFR: sx_buf_printf(b, "%d\t'Rtransf-%d_%d'\t%d\n", gt->index[i], st->index[gt->conts[i]], gt->paralog[i], gt->paralog[i]); break;
                case RGC: sx_buf_printf(b, "%d\t'Rgc-%d_%d'\t%d\n", gt->index[i], st->index[gt->conts[i]], gt->paralog[i], gt->paralog[i]); break;
                default:
                    if (names == NULL) sx_buf_printf(b, "%d\t'%d_%d'\t%d\n", gt->index[i], gt->sp_index[i], gt->paralog[i], gt->paralog[i]);
                    else sx_buf_printf(b, "%d\t'%s_%d'\t%d\n", gt->index[i], nm(names, gt->sp_index[i]), gt->paralog[i], gt->paralog[i]);
                    break;
                }
            }
        } else {
            switch (lt->kind[l]) {
            case LOSS: sx_buf_printf(b, "'Lost-%d_%d_%d'\t'Lost-%d_%d'\t%d\n", st->index[gt->conts[i]], gt->paralog[i], gt->replica[i], st->index[gt->conts[i]], gt->paralog[i], gt->paralog[i]); break;
            case RTRFR: sx_buf_printf(b, "'Rtransf-%d_%d_%d'\t'Rtransf-%d_%d'\t%d\n", st->index[gt->conts[i]], gt->paralog[i], gt->replica[i], st->index[gt->conts[i]], gt->paralog[i], gt->paralog[i]); break;
            case RGC: sx_buf_printf(b, "'Rgc-%d_%d_%d'\t'Rgc-%d_%d'\t%d\n", st->index[gt->conts[i]], gt->paralog[i], gt->replica[i], st->index[gt->conts[i]], gt->paralog[i], gt->paralog[i]); break;
            default:
                if (names == NULL) sx_buf_printf(b, "'%d_%d_%d'\t'%d_%d'\t%d\n", gt->sp_index[i], gt->paralog[i], gt->replica[i], gt->sp_index[i], gt->paralog[i], gt->paralog[i]);
                else sx_buf_printf(b, "'%s_%d_%d'\t'%s_%d'\t%d\n", nm(names, gt->sp_index[i]), gt->paralog[i], gt->replica[i], nm(names, gt->sp_index[i]), gt->paralog[i], gt->paralog[i]);
                break;
            }
        }
    }
    return SX_NO_ERROR;
}
