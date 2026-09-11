#include "sx_newick.h"
#include <ctype.h>

sx_names *sx_names_new(int n_names)
{
    sx_names *n = sx_calloc(1, sizeof(*n));
    n->n_names = n_names;
    n->names = sx_calloc((size_t)n_names + 1, sizeof(char *));
    n->names[0] = sx_malloc(sizeof("Internal node"));
    strcpy(n->names[0], "Internal node");
    return n;
}
void sx_names_free(sx_names *n)
{
    if (!n) return;
    for (int i = 0; i <= n->n_names; ++i) free(n->names[i]);
    free(n->names); free(n);
}

/* CheckNewickSTree port */
static long check_newick_stree(const char *tree)
{
    size_t i, k;
    int error = 0, n_left = 0, n_right = 0, n_branches = 0, n_nodes = 0, node = 0, f_leaf = 0, leaf = 0;
    long r_value = SX_NO_ERROR;
    char current = 0, previous, next;
    size_t len = strlen(tree);
    if (tree[0] != '(') error = 1;
    for (i = 0; i < len; ++i) {
        if (error == 1 || f_leaf == 1) {
            if (f_leaf == 1) { fprintf(stderr, "\nERROR: This internal node does not have more than one leaf, Is there a missed comma?\n"); f_leaf = 0; }
            else { fprintf(stderr, "\nERROR: There is something wrong in the tree\n"); error = 0; }
            for (k = 0; k + 1 <= i; k++) fprintf(stderr, "%c", tree[k]);
            fprintf(stderr, " <- HERE");
            r_value = SX_SETTINGS_ERROR;
        }
        current = tree[i];
        previous = tree[i > 0 ? i - 1 : 0];
        next = tree[i + 1];
        switch (current) {
        case '(':
            if (!(previous == '(' || previous == ',')) error = 1;
            if (next != '(') { ++n_nodes; leaf = 1; if (isalnum((unsigned char)next) == 0) error = 1; }
            ++n_left; node = 0; break;
        case ')':
            if (!(previous == ')' || isalnum((unsigned char)previous) != 0)) error = 1;
            if (!(next == ')' || next == ',' || next == ':' || next == '#' || next == '/' || next == '*' || next == '~')) error = 1;
            if (next == ';' && i == len - 2) error = 0;
            ++n_right; ++n_nodes;
            if (node == 0 && leaf == 1) f_leaf = 1;
            leaf = 0; node = 0; break;
        case ',':
            if (!(previous == ')' || isalnum((unsigned char)previous) != 0)) error = 1;
            if (next != '(') { ++n_nodes; leaf = 1; if (isalnum((unsigned char)next) == 0) error = 1; }
            node = 1; break;
        case ':': ++n_branches; /* fallthrough */
        case '*': case '#':
            if (!(previous == ')' || isalnum((unsigned char)previous) != 0)) error = 1;
            if (isdigit((unsigned char)next) == 0) error = 1;
            break;
        case '/':
            if (!(previous == ')' || isalnum((unsigned char)previous) != 0)) error = 1;
            if (isdigit((unsigned char)next) == 0) error = 1;
            if (leaf == 0) error = 1;
            break;
        default: break;
        }
    }
    if (current != ';') { fprintf(stderr, "\n\tTree is not finished with ;"); r_value = SX_SETTINGS_ERROR; }
    if (n_branches != n_nodes - 1) {
        fprintf(stderr, "\n\tThere are %s of nodes related with the number of branches. Maybe there are a missed comma, or the root node have branch length (it must not have it)\n", n_branches > n_nodes - 1 ? "a defect" : "an excess");
        r_value = SX_SETTINGS_ERROR;
    }
    if (n_left != n_right) { fprintf(stderr, "\n\tTree seems unbalanced (%d left and %d right parentheses)", n_left, n_right); r_value = SX_SETTINGS_ERROR; }
    return r_value;
}

static void refine_snodes(sx_stree *t, int node, int ind_persp, double gen_time)
{
    if (t->anc[node] != -1) {
        if (t->n_child[node] == 0) { if (t->n_replicas[node] == 0) t->n_replicas[node] = ind_persp; }
        else if (t->n_replicas[node] > 0) sx_check(SX_SETTINGS_ERROR, ": number of replicas per taxa applied to internal branches are not allowed\n");
        t->n_gen[node] = t->n_gen[t->anc[node]] + t->gen_length[node];
        t->time[node] = t->time[t->anc[node]] + t->gen_length[node] * t->gtime_mult[node] * gen_time;
    }
    for (int i = 0; i < t->n_child[node]; ++i) refine_snodes(t, SX_CHILD(t, node, i), ind_persp, gen_time);
}

static double check_ultrametricity(const sx_stree *t, int node)
{
    if (t->n_child[node] != 0) {
        double time = check_ultrametricity(t, SX_CHILD(t, node, 0));
        if (time == -1) return -1;
        for (int i = 1; i < t->n_child[node]; ++i) {
            double nt = check_ultrametricity(t, SX_CHILD(t, node, i));
            if (nt == -1 || fabs(time - nt) > FLOAT_PRECISION) return -1;
        }
        return time;
    }
    return t->time[node];
}

/* scan a run of characters satisfying pred into buf (bounded), advancing *pos; returns length */
static int scan_run(const char *s, size_t *pos, int (*pred)(char), char *buf, size_t cap)
{
    size_t n = 0;
    while (s[*pos] && pred(s[*pos]) && n + 1 < cap) { buf[n++] = s[*pos]; ++*pos; }
    buf[n] = 0;
    return (int)n;
}
static int is_digit_c(char c) { return c < 58 && c > 47; }
static int is_num_c(char c) { return c < 58 && (c > 47 || c == 46); }

sx_stree *sx_parse_newick_stree(const char *newick, sx_names **names_ptr, int verbosity, double gen_time, int Ne, double mu, int ind_persp)
{
    sx_check(check_newick_stree(newick), ": Error in the Newick species tree\n");
    if (gen_time == 0) sx_check(SX_UNEXPECTED_VALUE, ": Inproper generation time");
    size_t len = strlen(newick);
    char *buf = sx_malloc(len + 2);
    /* first pass: counts, replicas and maximum number of children */
    int n_leaves = 0, n_inodes = 0, n_gleaves = 0, n_priv = 0, max_children = 1;
    int depth = 0, *cnt = sx_calloc(len + 2, sizeof(int));
    for (size_t p = 0; p < len && newick[p] != ';'; ++p) {
        char c = newick[p];
        if ((c == '(' || c == ',') && newick[p + 1] != '(') ++n_leaves;
        if (c == ')') ++n_inodes;
        if (c == '(') { ++depth; cnt[depth] = 1; }
        else if (c == ',') { if (++cnt[depth] > max_children) max_children = cnt[depth]; }
        else if (c == ')') { --depth; }
        else if (c == '/') {
            size_t q = p + 1; int nr = 0;
            scan_run(newick, &q, is_digit_c, buf, len + 2);
            sscanf(buf, "%u", (unsigned *)&nr);
            n_gleaves += nr; ++n_priv;
            p = q - 1;
        }
    }
    free(cnt);
    if (max_children < 2) max_children = 2;
    n_gleaves += (n_leaves - n_priv) * ind_persp;
    sx_names *names = sx_names_new(n_leaves);
    *names_ptr = names;

    sx_stree *t = sx_stree_new(0, n_leaves, n_gleaves, max_children, gen_time, Ne, mu);
    int current = -1, anc = -1, root = -1, iteration = 0, ffree = 1;
    size_t step = 0;
    char code = newick[0];
    while (code != ';') {
        switch (code) {
        case '(':
            if (iteration != 0) {
                anc = current; current = sx_stree_add_node(t);
                if (t->n_child[anc] >= max_children) sx_fatal("Newick parser: too many children");
                SX_CHILD(t, anc, t->n_child[anc]) = current; ++t->n_child[anc]; t->anc[current] = anc;
            } else { current = sx_stree_add_node(t); root = current; }
            break;
        case ',': case ')':
            current = t->anc[current];
            anc = current >= 0 ? t->anc[current] : -1;
            break;
        case ':': { ++step; scan_run(newick, &step, is_num_c, buf, len + 2); --step; sscanf(buf, "%lf", &t->gen_length[current]); break; }
        case '*': { ++step; scan_run(newick, &step, is_num_c, buf, len + 2); --step; sscanf(buf, "%lf", &t->mu_mult[current]); break; }
        case '~': { ++step; scan_run(newick, &step, is_num_c, buf, len + 2); --step; sscanf(buf, "%lf", &t->gtime_mult[current]); break; }
        case '#': {
            ++step; scan_run(newick, &step, is_digit_c, buf, len + 2); --step;
            sscanf(buf, "%u", (unsigned *)&t->Ne_b[current]);
            if (t->Ne_b[current] == 0) sx_check(SX_SETTINGS_ERROR, ": Null branch specific population size is not allowed");
            break;
        }
        case '/': {
            ++step; scan_run(newick, &step, is_digit_c, buf, len + 2); --step;
            int nr = 0; sscanf(buf, "%u", (unsigned *)&nr);
            if (nr == 0) sx_check(SX_SETTINGS_ERROR, ": Null branch specific number of replicates is not allowed");
            t->n_replicas[current] = nr;
            break;
        }
        default: {
            size_t n = 0;
            while (code != '(' && code != ')' && code != ',' && code != ';' && code != ':' && code != '~' && n < (size_t)MAX_NAME) {
                buf[n++] = code; ++step; code = newick[step];
            }
            if (n == (size_t)MAX_NAME) fprintf(stderr, "Species names are being truncated, since they are bigger than the buffer. You could change this behaviour increasing the buffer size by changing the environment variable SIMPHY_MAXNAME\n");
            buf[n] = 0;
            --step;
            if (ffree > n_leaves) sx_fatal("Newick parser: more leaves than counted");
            names->names[ffree] = sx_malloc(n + 1); memcpy(names->names[ffree], buf, n + 1);
            anc = current; current = sx_stree_add_node(t);
            if (t->n_child[anc] >= max_children) sx_fatal("Newick parser: too many children");
            SX_CHILD(t, anc, t->n_child[anc]) = current; ++t->n_child[anc]; t->anc[current] = anc;
            t->sp_index[current] = ffree; ++ffree;
            break;
        }
        }
        ++step; code = newick[step]; ++iteration;
        if (iteration == MAX_IT) sx_check(SX_LOOP_ERROR, NULL);
    }
    free(buf);
    t->root = root;
    refine_snodes(t, root, ind_persp, gen_time);
    t->n_leaves = n_leaves; t->n_gleaves = n_gleaves;
    if (t->n_nodes != n_leaves + n_inodes) sx_fatal("Newick parser: node count mismatch (%d vs %d)", t->n_nodes, n_leaves + n_inodes);
    sx_stree_reindex(t, 1);
    if (check_ultrametricity(t, root) == -1) sx_check(SX_UNEXPECTED_VALUE, "The species tree is not ultrametric in time units. This doesn't allow to use the full model (transfers)\n");
    if (verbosity > 2) {
        sx_out("\n\t%d-node species tree correctly built", (n_leaves * 2) - 1);
        if (verbosity > 3) { sx_buf tb; sx_buf_init(&tb); sx_buf_puts(&tb, ": "); sx_write_stree(&tb, t, names, 0, 0); while (tb.len > 0 && tb.data[tb.len - 1] == '\n') --tb.len; tb.data[tb.len] = 0; sx_out("%s", tb.data); sx_buf_free(&tb); }
        sx_out("\n");
    }
    return t;
}
