#include "sx_newick.h"
#include <ctype.h>
#include "sx_ltree.h"

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

/* ======================= shared helpers ======================= */
static void print_xchar_error(const char *string, int x, const char *msg1, const char *msg)
{
    printf("%s", msg1);
    for (int i = 1; i < x; ++i) printf("%c", string[i]);
    printf(" %s \n", msg);
    fflush(stdout);
}

/* CheckNexusTree port: like CheckNewickSTree but allowing [&...] comments. */
static long check_nexus_tree(const char *tree)
{
    size_t i, k, len = strlen(tree);
    int error = 0, n_left = 0, n_right = 0, n_branches = 0, n_nodes = 0, node = 0, f_leaf = 0, leaf = 0, comment = 0;
    long r_value = SX_NO_ERROR;
    char current = 0, previous, next;
    if (tree[0] != '(') error = 1;
    for (i = 0; i < len; ++i) {
        if (error == 1 || f_leaf == 1) {
            if (f_leaf == 1) { fprintf(stderr, "\nERROR: This internal node does not have more than one leaf, Is there a missed comma?\n"); f_leaf = 0; }
            else { fprintf(stderr, "\nERROR: There is something wrong in the tree\n"); error = 0; }
            for (k = 0; k + 1 <= i; k++) fprintf(stderr, "%c", tree[k]);
            fprintf(stderr, " <- HERE");
            r_value = SX_SETTINGS_ERROR;
        }
        if (comment == 1) { if (tree[i] == ']') comment = 0; continue; }
        current = tree[i]; previous = tree[i > 0 ? i - 1 : 0]; next = tree[i + 1];
        switch (current) {
        case '(':
            if (!(previous == '(' || previous == ',')) error = 1;
            if (next != '(') { ++n_nodes; leaf = 1; if (isalnum((unsigned char)next) == 0) error = 1; }
            ++n_left; node = 0; break;
        case ')':
            if (!(previous == ')' || previous == ']' || isalnum((unsigned char)previous) != 0)) error = 1;
            if (!(next == ')' || next == ',' || next == ':' || next == '[')) error = 1;
            if (next == ';' && i == len - 2) error = 0;
            ++n_right; ++n_nodes;
            if (node == 0 && leaf == 1) f_leaf = 1;
            leaf = 0; node = 0; break;
        case ',':
            if (!(previous == ')' || previous == ']' || isalnum((unsigned char)previous) != 0)) error = 1;
            if (next != '(') { ++n_nodes; leaf = 1; if (isalnum((unsigned char)next) == 0) error = 1; }
            node = 1; break;
        case ':':
            ++n_branches;
            if (!(previous == ')' || previous == ']' || isalnum((unsigned char)previous) != 0)) error = 1;
            if (isdigit((unsigned char)next) == 0) error = 1;
            break;
        case '[':
            if (!(previous == ')' || isalnum((unsigned char)previous) != 0)) error = 1;
            if (next != '&') error = 1;
            comment = 1; break;
        default: break;
        }
    }
    if (current != ';') { fprintf(stderr, "\n\tTree is not finished with ;"); r_value = SX_SETTINGS_ERROR; }
    if (n_branches != n_nodes - 1) { fprintf(stderr, "\n\tThere are %s of nodes related with the number of branches. Maybe there are a missed comma, or the root node have branch length (it must not have it)\n", n_branches > n_nodes - 1 ? "a defect" : "an excess"); r_value = SX_SETTINGS_ERROR; }
    if (n_left != n_right) { fprintf(stderr, "\n\tTree seems unbalanced (%d left and %d right parentheses)", n_left, n_right); r_value = SX_SETTINGS_ERROR; }
    return r_value;
}

/* CheckNewickLTree port (locus flavour adds % and / codes). */
static long check_newick_ltree(const char *tree)
{
    size_t i, k, len = strlen(tree);
    int error = 0, n_left = 0, n_right = 0, n_branches = 0, n_nodes = 0, node = 0, f_leaf = 0, leaf = 0;
    long r_value = SX_NO_ERROR;
    char current = 0, previous, next;
    if (tree[0] != '(') error = 1;
    for (i = 0; i < len; ++i) {
        if (error == 1 || f_leaf == 1) {
            if (f_leaf == 1) { fprintf(stderr, "\nERROR: This internal node does not have more than one leaf, Is there a missed comma?\n"); f_leaf = 0; }
            else { fprintf(stderr, "\nERROR: There is something wrong in the tree\n"); error = 0; }
            for (k = 0; k + 1 <= i; k++) fprintf(stderr, "%c", tree[k]);
            fprintf(stderr, " <- HERE");
            r_value = SX_SETTINGS_ERROR;
        }
        current = tree[i]; previous = tree[i > 0 ? i - 1 : 0]; next = tree[i + 1];
        switch (current) {
        case '(':
            if (!(previous == '(' || previous == ',')) error = 1;
            if (next != '(') { ++n_nodes; leaf = 1; if (isalnum((unsigned char)next) == 0) error = 1; }
            ++n_left; node = 0; break;
        case ')':
            if (!(previous == ')' || isalnum((unsigned char)previous) != 0)) error = 1;
            if (!(next == ')' || next == ',' || next == ':' || next == '#' || next == '/' || next == '%' || next == '*' || next == '~')) error = 1;
            if (next == ';' && i == len - 2) error = 0;
            ++n_right; ++n_nodes;
            if (node == 0 && leaf == 1) f_leaf = 1;
            leaf = 0; node = 0; break;
        case ',':
            if (!(previous == ')' || isalnum((unsigned char)previous) != 0)) error = 1;
            if (next != '(') { ++n_nodes; leaf = 1; if (isalnum((unsigned char)next) == 0) error = 1; }
            node = 1; break;
        case ':': ++n_branches; /* fallthrough */
        case '#': case '*': case '%':
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
    if (n_branches != n_nodes - 1) { fprintf(stderr, "\n\tThere are %s of nodes related with the number of branches. Maybe there are a missed comma, or the root node have branch length (it must not have it)\n", n_branches > n_nodes - 1 ? "a defect" : "an excess"); r_value = SX_SETTINGS_ERROR; }
    if (n_left != n_right) { fprintf(stderr, "\n\tTree seems unbalanced (%d left and %d right parentheses)", n_left, n_right); r_value = SX_SETTINGS_ERROR; }
    return r_value;
}

/* parameters inside a Nexus comment; returns SX_NO_ERROR or an error code; *pos advanced past ']' position */
typedef struct { int *Ne, *kind, *paralog, *n_ind; double *mu_mult, *gtime_mult; } nexus_targets;
static long nexus_comment_params(const char *string, size_t *pos, nexus_targets *t, int locus)
{
    static const char *sp_params[4] = { "pop_size", "n_ind", "u_mult", "g_mult" };
    static const int sp_len[4] = { 8, 5, 6, 6 };
    static const char *lp_params[6] = { "pop_size", "kind_n", "paralog", "n_ind", "u_mult", "g_mult" };
    static const int lp_len[6] = { 8, 6, 7, 5, 6, 6 };
    const char **params = locus ? lp_params : sp_params;
    const int *plen = locus ? lp_len : sp_len;
    int n_params = locus ? 6 : 4;
    int i = 0, c_param = -1, out_param = 1;
    char buf[64];
    while (string[*pos] != ']' && string[*pos] != 0 && i <= MAX_IT) {
        if (out_param == 1) {
            c_param = -1;
            for (int k = 0; k < n_params; ++k) if (strncmp(params[k], string + *pos, (size_t)plen[k]) == 0) c_param = k;
            if (c_param == -1) { print_xchar_error(string, (int)*pos + 8, "\nNEXUS TREE PARSING ERROR\n", "|<- Unvalid param, please, revisit the manual to check the proper parameter spelling\n"); return SX_SETTINGS_ERROR; }
            else if (string[*pos + plen[c_param]] == '=') { out_param = 0; *pos += (size_t)plen[c_param] + 1; }
            else return SX_SETTINGS_ERROR;
        } else {
            size_t n = 0; char code = string[*pos];
            while (code < 58 && (code > 47 || code == 46) && n + 1 < sizeof buf) { buf[n++] = code; ++*pos; code = string[*pos]; }
            buf[n] = 0;
            int ok = 1;
            if (locus) {
                switch (c_param) {
                case 0: ok = sscanf(buf, "%d", t->Ne) == 1; break;
                case 1: ok = sscanf(buf, "%d", t->kind) == 1; break;
                case 2: ok = sscanf(buf, "%d", t->paralog) == 1; break;
                case 3: ok = sscanf(buf, "%d", t->n_ind) == 1; break;
                case 4: ok = sscanf(buf, "%lf", t->mu_mult) == 1; break;
                case 5: ok = sscanf(buf, "%lf", t->gtime_mult) == 1; break;
                default: return SX_SETTINGS_ERROR;
                }
            } else {
                switch (c_param) {
                case 0: ok = sscanf(buf, "%d", t->Ne) == 1; break;
                case 1: ok = sscanf(buf, "%d", t->n_ind) == 1; break;
                case 2: ok = sscanf(buf, "%lf", t->mu_mult) == 1; break;
                case 3: ok = sscanf(buf, "%lf", t->gtime_mult) == 1; break;
                default: return SX_SETTINGS_ERROR;
                }
            }
            if (!ok) { print_xchar_error(string, (int)*pos, "\nNEWICK PARSING ERROR\n", "|<- Unvalid value, please, check the NEXUS file\n"); return SX_SETTINGS_ERROR; }
            if (string[*pos] == ',') ++*pos;
            out_param = 1;
        }
        ++i;
    }
    if (i > MAX_IT) return SX_LOOP_ERROR;
    return SX_NO_ERROR;
}

/* ======================= species tree, Nexus flavour ======================= */
sx_stree *sx_parse_nexus_stree(const char *nexus, sx_names **names_ptr, int verbosity, double gen_time, int Ne, double mu, int ind_persp)
{
    sx_check(check_nexus_tree(nexus), ": Error in the Nexus species tree file\n");
    if (gen_time == 0) sx_check(SX_UNEXPECTED_VALUE, ": Inproper generation time");
    size_t len = strlen(nexus);
    char *buf = sx_malloc(len + 2);
    int n_leaves = 0, n_inodes = 0, n_gleaves = 0, max_children = 1, in_comment = 0, depth = 0;
    int *cnt = sx_calloc(len + 2, sizeof(int));
    for (size_t p = 0; p < len && nexus[p] != ';'; ++p) {
        char c = nexus[p];
        if (in_comment) { if (c == ']') in_comment = 0; continue; }
        if ((c == '(' || c == ',') && nexus[p + 1] != '(') ++n_leaves;
        if (c == ')') ++n_inodes;
        if (c == '[') { in_comment = 1; continue; }
        if (c == '(') { ++depth; cnt[depth] = 1; }
        else if (c == ',') { if (++cnt[depth] > max_children) max_children = cnt[depth]; }
        else if (c == ')') --depth;
    }
    free(cnt);
    if (max_children < 2) max_children = 2;
    sx_names *names = sx_names_new(n_leaves);
    *names_ptr = names;
    sx_stree *t = sx_stree_new(0, n_leaves, 0, max_children, gen_time, Ne, mu);
    int current = -1, anc = -1, root = -1, iteration = 0, ffree = 1;
    size_t step = 0;
    char code = nexus[0];
    while (code != ';') {
        switch (code) {
        case '(':
            if (iteration != 0) { anc = current; current = sx_stree_add_node(t); SX_CHILD(t, anc, t->n_child[anc]) = current; ++t->n_child[anc]; t->anc[current] = anc; }
            else { current = sx_stree_add_node(t); root = current; }
            break;
        case ',': case ')':
            current = t->anc[current]; anc = current >= 0 ? t->anc[current] : -1;
            break;
        case ':': { ++step; scan_run(nexus, &step, is_num_c, buf, len + 2); --step; sscanf(buf, "%lf", &t->gen_length[current]); break; }
        case '[': {
            step += 2;
            nexus_targets tg = { &t->Ne_b[current], NULL, NULL, &t->n_replicas[current], &t->mu_mult[current], &t->gtime_mult[current] };
            sx_check(nexus_comment_params(nexus, &step, &tg, 0), ": Error parsing Nexus tree comments as branch-specific parameters");
            if (t->n_replicas[current] > 1 && t->n_child[current] != 0) {
                print_xchar_error(nexus, (int)step, "\nNEWICK PARSING ERROR\n", "|<- Fixing the number of nodes/replicates in an internal node is not allowed, please, revisit the manual and your input trees\n");
                sx_check(SX_SETTINGS_ERROR, NULL);
            }
            if (t->n_replicas[current] > 1) n_gleaves += t->n_replicas[current] - ind_persp;
            break;
        }
        default: {
            size_t n = 0;
            while (code != '(' && code != ')' && code != ',' && code != ';' && code != ':' && code != '[' && n < (size_t)MAX_NAME) { buf[n++] = code; ++step; code = nexus[step]; }
            if (n == (size_t)MAX_NAME) fprintf(stderr, "Species names are being truncated, since they are bigger than the buffer. You could change this behaviour increasing the buffer size by changing the environment variable SIMPHY_MAXNAME\n");
            buf[n] = 0; --step;
            if (ffree > n_leaves) sx_fatal("Nexus parser: more leaves than counted");
            names->names[ffree] = sx_malloc(n + 1); memcpy(names->names[ffree], buf, n + 1);
            anc = current; current = sx_stree_add_node(t);
            n_gleaves += ind_persp;
            SX_CHILD(t, anc, t->n_child[anc]) = current; ++t->n_child[anc]; t->anc[current] = anc;
            t->sp_index[current] = ffree; ++ffree;
            break;
        }
        }
        ++step; code = nexus[step]; ++iteration;
        if (iteration == MAX_IT) sx_check(SX_LOOP_ERROR, NULL);
    }
    free(buf);
    t->root = root;
    refine_snodes(t, root, ind_persp, gen_time);
    t->n_leaves = n_leaves; t->n_gleaves = n_gleaves;
    if (t->n_nodes != n_leaves + n_inodes) sx_fatal("Nexus parser: node count mismatch");
    sx_stree_reindex(t, 1);
    if (check_ultrametricity(t, root) == -1) sx_check(SX_UNEXPECTED_VALUE, "The species tree is not ultrametric in time units. This doesn't allow to use the full model (transfers)\n");
    if (verbosity > 2) {
        sx_out("\n\t%d-node species tree correctly built", (n_leaves * 2) - 1);
        if (verbosity > 3) { sx_buf tb; sx_buf_init(&tb); sx_buf_puts(&tb, ": "); sx_write_stree(&tb, t, names, 0, 0); while (tb.len > 0 && tb.data[tb.len - 1] == '\n') --tb.len; tb.data[tb.len] = 0; sx_out("%s", tb.data); sx_buf_free(&tb); }
        sx_out("\n");
    }
    return t;
}

/* ======================= locus trees ======================= */
static void refine_lnodes(sx_ltree *t, int node, int ind_persp, double gen_time, int *n_dup, int *n_loss, int *n_trans, int *n_gc)
{
    switch (t->kind[node]) {
    case DUP: (*n_dup)++; break;
    case LOSS: (*n_loss)++; __attribute__((fallthrough));
    case GC: (*n_gc)++; break;
    case TRFR: (*n_trans)++; break;
    default: break;
    }
    if (t->anc[node] != -1) {
        if (t->n_child[node] == 0) {
            switch (t->kind[node]) {
            case LOSS: case RTRFR: case RGC: t->n_lin[node] = 0; break;
            default: if (t->n_lin[node] == 0) t->n_lin[node] = ind_persp; break;
            }
            t->n_ilin[node] = t->n_lin[node];
        }
        t->n_gen[node] = t->n_gen[t->anc[node]] + t->gen_length[node];
        t->time[node] = t->time[t->anc[node]] + t->gen_length[node] * t->gtime_mult[node] * gen_time;
    }
    for (int i = 0; i < t->n_child[node]; ++i) refine_lnodes(t, SX_CHILD(t, node, i), ind_persp, gen_time, n_dup, n_loss, n_trans, n_gc);
}

static double check_ultrametricity_l(const sx_ltree *t, int node)
{
    if (t->n_child[node] != 0) {
        double time = check_ultrametricity_l(t, SX_CHILD(t, node, 0));
        if (time == -1) return -1;
        for (int i = 1; i < t->n_child[node]; ++i) {
            int c = SX_CHILD(t, node, i);
            if (t->kind[c] != LOSS && t->kind[c] != RTRFR && t->kind[c] != RGC) {
                double nt = check_ultrametricity_l(t, c);
                if (nt == -1 || fabs(time - nt) > FLOAT_PRECISION) return -1;
            }
        }
        return time;
    }
    return t->time[node];
}

/* shared body: nexus=0 → extended Newick codes, nexus=1 → [&...] comments */
static sx_ltree *parse_ltree(const char *str, int nexus, sx_names **names_ptr, int verbosity, double gen_time, int Ne, double mu, int ind_persp, int *n_dup, int *n_loss, int *n_trans, int *n_gc)
{
    sx_check(nexus ? check_nexus_tree(str) : check_newick_ltree(str), nexus ? ": Nexus locus tree error" : ": Newick locus tree error");
    if (gen_time == 0) sx_check(SX_UNEXPECTED_VALUE, ": Inproper generation time");
    size_t len = strlen(str);
    char *buf = sx_malloc(len + 2);
    int n_leaves = 0, n_inodes = 0, n_gleaves = 0, n_priv = 0, max_children = 1, in_comment = 0, depth = 0;
    int *cnt = sx_calloc(len + 2, sizeof(int));
    for (size_t p = 0; p < len && str[p] != ';'; ++p) {
        char c = str[p];
        if (in_comment) { if (c == ']') in_comment = 0; continue; }
        if ((c == '(' || c == ',') && str[p + 1] != '(') ++n_leaves;
        if (c == ')') ++n_inodes;
        if (nexus && c == '[') { in_comment = 1; continue; }
        if (c == '(') { ++depth; cnt[depth] = 1; }
        else if (c == ',') { if (++cnt[depth] > max_children) max_children = cnt[depth]; }
        else if (c == ')') --depth;
        else if (!nexus && c == '/') { size_t q = p + 1; int nr = 0; scan_run(str, &q, is_digit_c, buf, len + 2); sscanf(buf, "%u", (unsigned *)&nr); n_gleaves += nr; ++n_priv; p = q - 1; }
    }
    free(cnt);
    if (max_children < 2) max_children = 2;
    int n_tleaves = n_leaves;
    sx_names *names = sx_names_new(n_leaves);
    *names_ptr = names;
    sx_ltree *t = sx_ltree_new(0, n_leaves, 0, max_children, gen_time, Ne, mu);
    int current = -1, anc = -1, root = -1, iteration = 0, ffree = 1;
    size_t step = 0;
    char code = str[0];
    while (code != ';') {
        switch (code) {
        case '(':
            if (iteration != 0) { anc = current; current = sx_ltree_add_node(t); SX_CHILD(t, anc, t->n_child[anc]) = current; ++t->n_child[anc]; t->anc[current] = anc; }
            else { current = sx_ltree_add_node(t); root = current; }
            break;
        case ',': case ')':
            current = t->anc[current]; anc = current >= 0 ? t->anc[current] : -1;
            break;
        case ':': { ++step; scan_run(str, &step, is_num_c, buf, len + 2); --step; sscanf(buf, "%lf", &t->gen_length[current]); break; }
        case '*': if (nexus) goto name; { ++step; scan_run(str, &step, is_num_c, buf, len + 2); --step; sscanf(buf, "%lf", &t->mu_mult[current]); break; }
        case '~': if (nexus) goto name; { ++step; scan_run(str, &step, is_num_c, buf, len + 2); --step; sscanf(buf, "%lf", &t->gtime_mult[current]); break; }
        case '#': if (nexus) goto name; { ++step; scan_run(str, &step, is_digit_c, buf, len + 2); --step; sscanf(buf, "%u", (unsigned *)&t->Ne_b[current]); if (t->Ne_b[current] == 0) sx_check(SX_SETTINGS_ERROR, ": Null branch specific population size is not allowed"); break; }
        case '/': if (nexus) goto name; { ++step; scan_run(str, &step, is_digit_c, buf, len + 2); --step; int nr = 0; sscanf(buf, "%u", (unsigned *)&nr); if (nr == 0) sx_check(SX_SETTINGS_ERROR, ": Null number of replicates is not allowed"); t->n_lin[current] = nr; break; }
        case '%': if (nexus) goto name; { ++step; scan_run(str, &step, is_digit_c, buf, len + 2); --step; int kd = 0; sscanf(buf, "%u", (unsigned *)&kd); t->kind[current] = kd; if (kd == LOSS || kd == RGC || kd == RTRFR) --n_tleaves; break; }
        case '_': if (nexus) goto name; { ++step; scan_run(str, &step, is_digit_c, buf, len + 2); --step; int pl = 0; sscanf(buf, "%u", (unsigned *)&pl); t->paralog[current] = pl; break; }
        case '[': {
            if (!nexus) goto name;
            step += 2;
            nexus_targets tg = { &t->Ne_b[current], &t->kind[current], &t->paralog[current], &t->n_lin[current], &t->mu_mult[current], &t->gtime_mult[current] };
            int before = t->n_lin[current];
            sx_check(nexus_comment_params(str, &step, &tg, 1), ": parsing Nexus tree comments");
            if (t->n_lin[current] != before) t->n_ilin[current] = t->n_lin[current];
            if (t->n_lin[current] > 0 && t->n_lin[current] != ind_persp) {
                n_gleaves += t->n_lin[current]; ++n_priv;
                if (t->n_child[current] != 0) { print_xchar_error(str, (int)step, "\nNEWICK PARSING ERROR\n", "|<- Fixing the number of nodes/replicates in an internal node is not allowed, please, revisit the manual and your input trees\n"); sx_check(SX_SETTINGS_ERROR, NULL); }
            }
            if (t->kind[current] == LOSS || t->kind[current] == RGC || t->kind[current] == RTRFR) --n_tleaves;
            break;
        }
        default:
        name: {
            size_t n = 0;
            if (nexus) { while (code != '(' && code != ')' && code != ',' && code != ';' && code != ':' && code != '[' && n < (size_t)MAX_NAME) { buf[n++] = code; ++step; code = str[step]; } }
            else { while (code != '(' && code != ')' && code != ',' && code != ';' && code != ':' && code != '_' && code != '~' && n < (size_t)MAX_NAME) { buf[n++] = code; ++step; code = str[step]; } }
            if (n == (size_t)MAX_NAME) fprintf(stderr, "%sSpecies names are being truncated, since they are bigger than the buffer. You could change this behaviour increasing the buffer size by changing the environment variable SIMPHY_MAXNAME\n", nexus ? "\n" : "");
            buf[n] = 0; --step;
            if (ffree > n_leaves) sx_fatal("locus tree parser: more leaves than counted");
            names->names[ffree] = sx_malloc(n + 1); memcpy(names->names[ffree], buf, n + 1);
            anc = current; current = sx_ltree_add_node(t);
            SX_CHILD(t, anc, t->n_child[anc]) = current; ++t->n_child[anc]; t->anc[current] = anc;
            t->sp_index[current] = ffree; ++ffree;
            break;
        }
        }
        ++step; code = str[step]; ++iteration;
        if (iteration == MAX_IT) sx_check(SX_LOOP_ERROR, NULL);
    }
    free(buf);
    n_gleaves += (n_tleaves - n_priv) * ind_persp;
    t->root = root;
    refine_lnodes(t, root, ind_persp, gen_time, n_dup, n_loss, n_trans, n_gc);
    t->n_leaves = n_leaves; t->n_gleaves = n_gleaves;
    if (t->n_nodes != n_leaves + n_inodes) sx_fatal("locus tree parser: node count mismatch");
    t->species_tree = NULL;
    /* PostReorderLNodes: post-order indices before the echo (the caller's collapse keeps them) */
    {
        int *order = sx_malloc(sizeof(int) * (size_t)t->n_nodes);
        int n = sx_order_traversal(t->n_child, t->child, t->max_children, root, 1, order, t->n_nodes);
        for (int k = 0; k < n; ++k) t->index[order[k]] = k;
        free(order);
    }
    if (check_ultrametricity_l(t, root) == -1) sx_check(SX_UNEXPECTED_VALUE, nexus ? ": locus tree is not ultrametric in time units. This doesn't allow to use the full model (transfers)\n" : "The species tree is not ultrametric in time units. This doesn't allow to use the full model (transfers)\n");
    if (verbosity > 4) sx_out("\n\t\tCalculating lineage count probabilities of entering and leaving each locus tree branch... ");
    sx_calc_probs_nlineages(t, root, Ne, 0, verbosity);
    if (verbosity > 4) sx_out("\n\t\tDone\n");
    if (verbosity > 2) {
        sx_out(nexus ? "\n\t%d-node locus tree correctly built" : "\n\t\t %d-node locus tree correctly built", (n_leaves * 2) - 1);
        if (verbosity > 3) { sx_buf tb; sx_buf_init(&tb); sx_buf_puts(&tb, ": "); sx_write_ltree(&tb, t, names, 0, 0); while (tb.len > 0 && tb.data[tb.len - 1] == '\n') --tb.len; tb.data[tb.len] = 0; sx_out("%s", tb.data); sx_buf_free(&tb); }
        sx_out("\n");
    }
    /* CollapseLTree(post-order, probs kept) */
    sx_ltree_collapse(t, t->n_nodes, n_leaves, 1);
    return t;
}

sx_ltree *sx_parse_newick_ltree(const char *newick, sx_names **names, int verbosity, double gen_time, int Ne, double mu, int ind_persp, int *n_dup, int *n_loss, int *n_trans, int *n_gc)
{ return parse_ltree(newick, 0, names, verbosity, gen_time, Ne, mu, ind_persp, n_dup, n_loss, n_trans, n_gc); }
sx_ltree *sx_parse_nexus_ltree(const char *nexus, sx_names **names, int verbosity, double gen_time, int Ne, double mu, int ind_persp, int *n_dup, int *n_loss, int *n_trans, int *n_gc)
{ return parse_ltree(nexus, 1, names, verbosity, gen_time, Ne, mu, ind_persp, n_dup, n_loss, n_trans, n_gc); }

/* ======================= Nexus files ======================= */
static int first_noblank(const char *s) { int i = 0; while (s[i] && isblank((unsigned char)s[i])) ++i; return i; }

long sx_nexus_count(FILE *f, int *n_trees)
{
    char *line = NULL; size_t cap = 0; ssize_t n;
    int in_trees = 0;
    *n_trees = 0;
    rewind(f);
    if ((n = getline(&line, &cap, f)) == -1 || strcmp(line, "#NEXUS\n") != 0) { fprintf(stderr, "Improper Nexus trees input file\n"); free(line); return SX_SETTINGS_ERROR; }
    while ((n = getline(&line, &cap, f)) != -1) {
        if (in_trees == 1) {
            char *nl = strchr(line, '\n'); if (nl) *nl = 0;
            if (strcmp(line, "End;") == 0 || strcmp(line, "END;") == 0 || strcmp(line, "end;") == 0) break;
            char *p = line + first_noblank(line);
            char *sp = strchr(p, ' '); if (sp) *sp = 0;
            if (strcmp(p, "tree") == 0) ++*n_trees;
            else { fprintf(stderr, "Improper Nexus trees input file: Error in line %s\n", line); free(line); return SX_IO_ERROR; }
        } else {
            char up[17]; int j;
            for (j = 0; j < 16 && line[j]; ++j) up[j] = (char)toupper((unsigned char)line[j]);
            up[j] = 0;
            if (strcmp(up, "BEGIN TREES;\n") == 0) in_trees = 1;
        }
    }
    free(line);
    return SX_NO_ERROR;
}

long sx_nexus_init(FILE *f)
{
    char *line = NULL; size_t cap = 0; ssize_t n;
    rewind(f);
    while ((n = getline(&line, &cap, f)) != -1) {
        for (ssize_t j = 0; j < n; ++j) line[j] = (char)toupper((unsigned char)line[j]);
        if (strcmp(line, "BEGIN TREES;\n") == 0) { free(line); return SX_NO_ERROR; }
    }
    free(line);
    return SX_IO_ERROR;
}

long sx_nexus_next(FILE *f, char **tree_str)
{
    char *line = NULL; size_t cap = 0;
    if (getline(&line, &cap, f) == -1) { fprintf(stderr, "Error getting next NEXUS tree, line: %s\n", line ? line : ""); free(line); return SX_IO_ERROR; }
    char *p = strchr(line, '(');
    if (!p) { fprintf(stderr, "Error parsing Nexus tree, no parenthesis present!!!\n"); free(line); return SX_LOOP_ERROR; }
    char *nl = strchr(p, '\n'); if (nl) *nl = 0;
    free(*tree_str);
    *tree_str = sx_malloc(strlen(p) + 1); strcpy(*tree_str, p);
    free(line);
    return SX_NO_ERROR;
}
