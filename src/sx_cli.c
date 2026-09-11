#include "sx_cli.h"
#include <ctype.h>
#include <time.h>

sx_settings *sx_settings_new(void)
{
    sx_settings *s = sx_calloc(1, sizeof(*s));
    s->ns_trees = 1; s->ng_trees = 1; s->verbosity = 1; s->min_lleaves = 1; s->min_lsleaves = 0; s->stats = 0; s->map = 0; s->db = 0; s->params = 0;
    s->command = 1; s->t_kind = 1; s->out_time = 0; s->out_daughters = 0; s->out_labels = 0;
    s->epsilon_brent = 0.000001f;
    s->u_seed = (unsigned long)(time(NULL) * clock());
    s->mode = SX_MODE_COMPAT; s->threads = 0;
    sx_set_double(&s->alpha_s, 0);
    sx_set_uint(&s->Ne, 0); sx_set_double(&s->mu, 1); sx_set_double(&s->gen_time, 1); sx_set_uint(&s->nl_trees, 1);
    sx_set_double(&s->sb_rate, 0); sx_set_double(&s->sd_rate, 0); sx_set_uint(&s->ind_per_sp, 1); sx_set_double(&s->outgroup, 0);
    sx_set_double(&s->bds_length, 0); sx_set_uint(&s->bds_leaves, 0);
    sx_set_double(&s->b_rate, 0); sx_set_double(&s->d_rate, 0); sx_set_double(&s->t_rate, 0); sx_set_double(&s->gc_rate, 0); sx_set_double(&s->alpha_l, 0);
    sx_set_double(&s->salpha_g, 0);
    sx_set_pointer_double(&s->lb_rate, &s->b_rate); sx_set_pointer_double(&s->ld_rate, &s->d_rate);
    sx_set_pointer_double(&s->lt_rate, &s->t_rate); sx_set_pointer_double(&s->lgc_rate, &s->gc_rate);
    sx_set_pointer_double(&s->lalpha_g, &s->salpha_g);
    sx_set_pointer_double(&s->alpha_g, &s->lalpha_g);
    const char *names[23] = { "SP", "SL", "SI", "RL", "GB", "GD", "GT", "GG", "HH", "GP", "LB", "LD", "LT", "LG", "SB", "SD", "ST", "SO", "SG", "SU", "HS", "HL", "HG" };
    sx_sunit *ptrs[23] = { &s->Ne, &s->bds_leaves, &s->ind_per_sp, &s->nl_trees, &s->b_rate, &s->d_rate, &s->t_rate, &s->gc_rate, &s->lalpha_g, &s->salpha_g, &s->lb_rate, &s->ld_rate, &s->lt_rate, &s->lgc_rate, &s->sb_rate, &s->sd_rate, &s->bds_length, &s->outgroup, &s->gen_time, &s->mu, &s->alpha_s, &s->alpha_l, &s->alpha_g };
    for (int i = 0; i < 23; ++i) { strcpy(s->table_storage[i].name, names[i]); s->table_storage[i].p = ptrs[i]; }
    s->table.table = s->table_storage; s->table.n_duples = 23;
    return s;
}

void sx_settings_free(sx_settings *s)
{
    if (!s) return;
    free(s->species_tree_str); free(s->locus_tree_str); free(s->out_name); free(s->stree_iname); free(s->ltree_iname); free(s->confile_text);
    if (s->stree_ifile) fclose(s->stree_ifile);
    if (s->ltree_ifile) fclose(s->ltree_ifile);
    free(s);
}

void sx_print_usage(void)
{
    printf("\nUsage: simphy-x -[Parameter code] value(i|r|c|b|*) ...\n\nValue kinds\n\ti=integer\n\tr=real\n\tc=character string\n\tb=boolean(1 or 0)\n\t*=sampling notation\n\nPARAMETERS\n__________\n__________\n\n-R...→ Replicates\n_________________\n\n\t-RS i: Number of species tree replicates (study replicates).\n\t-RL *: Number of locus trees per species tree.\n\t-RG i: Number of gene trees per locus tree (Not for general usage).\n\n-G... → Genome-wide parameters (sampled for each species tree)\n______________________________________________________________\n\n\t-GB *: Duplication parameter (to use with LB)\n\t-GD *: Loss rate parameter (to use with LD).\n\t-GT *: Transfer parameter (to use with LT).\n\t-GG *: Gene conversion parameter (to use with LG).\n\t-GP *: Gene-by-lineage-specific parameter (to use with HG).\n\n-S... → Species tree\n____________________\n\n\t-S c: Fixed species tree (extended Newick format).\n\t-SR c: Nexus file with species trees.\n\t-SB *: Speciation rate (events/time unit).\n\t-SD *: Extinction rate (events/time unit).\n\t-ST *: Species tree height (time units).\n\t-SL *: Number of taxa.\n\t-SO *: Ratio between ingroup height and the branch from the root to the ingroup. If this parameter is not set the outgroup is not simulated.\n\t-SI *: Number of individuals per species.\n\t-SP *: Tree-wide effective population size.\n\t-SU *: Tree-wide substitution rate.\n\t-SG *: Tree-wide generation time.\n\n-L... → Locus tree\n__________________\n\n\t-L c: Locus tree (extended Newick format).\n\t-LR c: Nexus file containing locus trees.\n\t-LB *: Duplication rate (events/generation).\n\t-LD *: Loss rate (events/generation).\n\t-LT *: Horizontal gene transfer (HGT) rate (events/generation). \n\t-LG *: Gene conversion (GC) rate (events/generation).\n\t-LK b: Distance-dependent HGT/GC: Determines whether the sampling of receptors of genetic material depends (1) or not (0) on the evolutionary distance between candidates and donors.\n\t-LL i: Minimum number of locus tree leaves.\n\t-LS i: Minimum number of species represented by the locus tree leaves.\n\n-H...→ Substitution rate heterogeneity parameters\n_________________________________________________\n\n\t-HS *: Species-specific branch rate heterogeneity modifiers.\n\t-HL *: Gene-family-specific rate heterogeneity modifiers.\n\t-HH *:Gene-by-lineage-specific locus tree parameter (to use with the HG argument below).\n\t-HG *: Gene-by-lineage-specific rate heterogeneity modifiers.\n\n-C... → Global options\n______________________\n\n\t-CS i: Random number generator seed.\n\t-CE r: Precision of the Brent’s method for root-finding when sampling the multilocus coalescent. (Not for general usage)\n\n-I c: Input configuration file\n\n-V [0,6]: Verbosity. (Note: the bigger the verbosity the slower SimPhy becomes. Levels over 3 may only make sense for debugging or study the way SimPhy works).\n\n-O...→ Output\n_____________\n\n\t-O c: Common output prefix-name (for folder and names).\n\t-OT b: Determines whether the species and locus tree branches are written in number of generations (0) or time units (1).\n\t-OM b: Activates the tree mapping output.\n\t-OD b: Activates the SQLite database output.\n\t-OP b: Activates the logging of sampled options.\n\t-OC b: Activates the logging of original command line parameters and input configuration files.\n\t-OL b: Activates the output of trees with internal nodes labelled by its post-order id starting from 0.\n\t-ON b: Activates the output of the bounded locus subtrees file.\n\n-X...→ simphy-x options\n_______________________\n\n\t-XM c: Mode: compat (default; byte-identical to SimPhy 1.0.2) or fast (same model, independent per-locus random streams, no arbitrary precision).\n\t-XT i: Number of threads (default: all available cores).\n\t-XV: Print the version and exit.\n\n\nSampling notation\n_________________\n\nNotation squeme= Distribution_code:parameter_1,parameter_2,...,parameter_n\nDistribution codes:\n\tF: fixed value\n\tU: Uniform\n\tN: Normal\n\tE: Exponential\n\tG: Gamma\n\tLN: Lognormal\n\tLU: LogUniform\n\tSL: Lognormal * constant\nExample: N:1,1 (Normal with mean and sd equal 1)\n\nExample\n_______\n\nsimphy-x -sb f:0.000001 -ld f:0.0000005 -lb f:0.0000005 -lt f:0.0000005 -rs 100 -rl U:10,100 -rg 1 -o SimPhy_test -sp f:10000 -su f:0.00001 -sg f:1 -sl U:20,50 -st f:1000000 -om 1 -v 2 -od 1 -op 1 -oc 1 -on 1 -cs 22\n\n");
    fflush(stdout);
}

static void print_xstring_error(char **argv, int x, const char *msg)
{
    printf("\n\nSETTINGS ERROR\n");
    for (int i = 1; i < x; ++i) printf("%s ", argv[i]);
    printf(" %s \n", msg);
    fflush(stdout);
}

#define ARG_ERR() do { print_xstring_error(argv, i, "|<– ERROR in this parameter\n"); return SX_SETTINGS_ERROR; } while (0)
#define ARG_UNK() do { print_xstring_error(argv, i, "|<– Unrecognized parameter\n"); return SX_SETTINGS_ERROR; } while (0)
#define PARSE_SU(dest) do { if (sx_parse_sampling(val, &(dest), s->table) != SX_NO_ERROR) ARG_ERR(); } while (0)
#define PARSE_UINT(dest) do { if (sscanf(val, "%u", (unsigned *)&(dest)) == 0) ARG_ERR(); } while (0)

/* Apply one option (code letters after '-', value string). argv/i are used for error messages. */
static long apply_option(sx_settings *s, const char *arg, const char *val, char **argv, int i)
{
    char code = arg[1], sub_code = arg[2];
    switch (toupper((unsigned char)code)) {
    case 'R':
        switch (toupper((unsigned char)sub_code)) {
        case 'S': PARSE_UINT(s->ns_trees); if (s->ns_trees < 1) ARG_ERR(); break;
        case 'G': PARSE_UINT(s->ng_trees); if (s->ng_trees < 1) ARG_ERR(); break;
        case 'L': PARSE_SU(s->nl_trees); break;
        default: ARG_UNK();
        }
        break;
    case 'H':
        switch (toupper((unsigned char)sub_code)) {
        case 'S': PARSE_SU(s->alpha_s); break;
        case 'L': PARSE_SU(s->alpha_l); break;
        case 'G': PARSE_SU(s->alpha_g); break;
        case 'H': PARSE_SU(s->lalpha_g); break;
        default: ARG_UNK();
        }
        break;
    case 'V':
        PARSE_UINT(s->verbosity);
        break;
    case 'O':
        switch (toupper((unsigned char)sub_code)) {
        case '\0': free(s->out_name); s->out_name = sx_calloc(strlen(val) + 1, 1); sscanf(val, "%s", s->out_name); break;
        case 'T': PARSE_UINT(s->out_time); break;
        case 'S': PARSE_UINT(s->stats); break;
        case 'M': PARSE_UINT(s->map); break;
        case 'D': PARSE_UINT(s->db); break;
        case 'P': PARSE_UINT(s->params); break;
        case 'L': PARSE_UINT(s->out_labels); break;
        case 'N': PARSE_UINT(s->out_daughters); break;
        case 'C': PARSE_UINT(s->command); break;
        default: ARG_UNK();
        }
        break;
    case 'S':
        switch (toupper((unsigned char)sub_code)) {
        case '\0': free(s->species_tree_str); s->species_tree_str = sx_calloc(strlen(val) + 1, 1); sscanf(val, "%s", s->species_tree_str); break;
        case 'R':
            print_xstring_error(argv, i, "|<- Nexus species tree files (-SR) are not supported yet by simphy-x\n");
            return SX_SETTINGS_ERROR;
        case 'B': PARSE_SU(s->sb_rate); break;
        case 'D': PARSE_SU(s->sd_rate); break;
        case 'T': PARSE_SU(s->bds_length); break;
        case 'L': PARSE_SU(s->bds_leaves); break;
        case 'O': PARSE_SU(s->outgroup); break;
        case 'I': PARSE_SU(s->ind_per_sp); break;
        case 'G': PARSE_SU(s->gen_time); break;
        case 'P': PARSE_SU(s->Ne); break;
        case 'U': PARSE_SU(s->mu); break;
        default: ARG_UNK();
        }
        break;
    case 'L':
        switch (toupper((unsigned char)sub_code)) {
        case '\0': free(s->locus_tree_str); s->locus_tree_str = sx_calloc(strlen(val) + 1, 1); sscanf(val, "%s", s->locus_tree_str); break;
        case 'R':
            print_xstring_error(argv, i, "|<- Nexus locus tree files (-LR) are not supported yet by simphy-x\n");
            return SX_SETTINGS_ERROR;
        case 'B': PARSE_SU(s->lb_rate); break;
        case 'D': PARSE_SU(s->ld_rate); break;
        case 'T': PARSE_SU(s->lt_rate); break;
        case 'G': PARSE_SU(s->lgc_rate); break;
        case 'K': PARSE_UINT(s->t_kind); if (s->t_kind < 0 || s->t_kind > 1) ARG_ERR(); break;
        case 'L': PARSE_UINT(s->min_lleaves); if (s->min_lleaves < 2) ARG_ERR(); break;
        case 'S': PARSE_UINT(s->min_lsleaves); if (s->min_lsleaves < 2) ARG_ERR(); break;
        default: ARG_UNK();
        }
        break;
    case 'G':
        switch (toupper((unsigned char)sub_code)) {
        case 'B': PARSE_SU(s->b_rate); break;
        case 'D': PARSE_SU(s->d_rate); break;
        case 'T': PARSE_SU(s->t_rate); break;
        case 'G': PARSE_SU(s->gc_rate); break;
        case 'P': PARSE_SU(s->salpha_g); break;
        default: ARG_UNK();
        }
        break;
    case 'C':
        switch (toupper((unsigned char)sub_code)) {
        case 'S': if (sscanf(val, "%lu", &s->u_seed) == 0) ARG_ERR(); break;
        case 'E': if (sscanf(val, "%f", &s->epsilon_brent) == 0 || s->epsilon_brent < 0) ARG_ERR(); break;
        default: ARG_UNK();
        }
        break;
    case 'X':
        switch (toupper((unsigned char)sub_code)) {
        case 'M':
            if (strcasecmp(val, "compat") == 0) s->mode = SX_MODE_COMPAT;
            else if (strcasecmp(val, "fast") == 0) s->mode = SX_MODE_FAST;
            else ARG_ERR();
            break;
        case 'T': PARSE_UINT(s->threads); break;
        default: ARG_UNK();
        }
        break;
    default:
        ARG_UNK();
    }
    return SX_NO_ERROR;
}

/* Configuration file: one "-CODE value" per line; '#' and "//" start comments. */
static long apply_config_file(sx_settings *s, const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) { perror("Error opening input_file: "); return SX_IO_ERROR; }
    char *line = NULL; size_t cap = 0; ssize_t n;
    int lineno = 0;
    while ((n = getline(&line, &cap, f)) != -1) {
        ++lineno;
        char *p = line;
        while (*p == ' ' || *p == '\t') ++p;
        if (*p == '#' || (p[0] == '/' && p[1] == '/') || *p == '\n' || *p == '\0') continue;
        if (*p != '-') { fprintf(stderr, "Error in the parameter: %s\n", line); free(line); fclose(f); return SX_SETTINGS_ERROR; }
        char opt[8] = { 0 }; int off = 0;
        if (sscanf(p, "%7s%n", opt, &off) != 1) { fprintf(stderr, "Error in the parameter: %s\n", line); free(line); fclose(f); return SX_SETTINGS_ERROR; }
        char *val = p + off;
        while (*val == ' ' || *val == '\t') ++val;
        char *e = val + strlen(val);
        while (e > val && (e[-1] == '\n' || e[-1] == '\r' || e[-1] == ' ' || e[-1] == '\t')) *--e = 0;
        if (toupper((unsigned char)opt[1]) == 'I' && opt[2] == 0) continue;   /* nested -I ignored */
        char *fake[2] = { (char *)path, opt };
        long rc = apply_option(s, opt, val, fake, 1);
        if (rc != SX_NO_ERROR) { free(line); fclose(f); return rc; }
    }
    free(line); fclose(f);
    return SX_NO_ERROR;
}

long sx_get_settings(sx_settings *s, int argc, char **argv)
{
    int i;
    s->argc = argc; s->argv = argv;
    if (argc < 2) { sx_print_usage(); return SX_TERMINATE_NOERROR; }
    /* configuration file first: command-line options take precedence (as in SimPhy) */
    for (i = 1; i + 1 < argc; ++i) {
        if (argv[i][0] == '-' && toupper((unsigned char)argv[i][1]) == 'I' && argv[i][2] == '\0') {
            if (argc > 3) printf("\n\tWARNING! You are using both the input file option -I and command line parameters. If you define a parameter using both input options, command line options will have preference\n");
            s->confile_name = argv[i + 1];
            {
                FILE *cf = fopen(argv[i + 1], "r");
                if (cf) {
                    size_t cap = 4096, len = 0; char *txt = sx_malloc(cap); int c;
                    while ((c = fgetc(cf)) != EOF) { if (len + 2 > cap) { cap *= 2; txt = sx_realloc(txt, cap); } txt[len++] = (char)c; }
                    txt[len] = 0; fclose(cf); free(s->confile_text); s->confile_text = txt;
                }
            }
            long rc = apply_config_file(s, argv[i + 1]);
            if (rc != SX_NO_ERROR) { print_xstring_error(argv, i + 1, "|<- ERROR in this parameter\n"); return rc; }
        }
    }
    for (i = 1; i < argc; ++i) {
        const char *arg = argv[i];
        if (arg[0] != '-') { print_xstring_error(argv, i, "|<– ERROR in this parameter\n"); return SX_SETTINGS_ERROR; }
        if (arg[1] == '-') {
            if (strcmp(arg, "--help") == 0) { sx_print_usage(); return SX_TERMINATE_NOERROR; }
            if (strcmp(arg, "--version") == 0) { printf("simphy-x %s (SimPhy 1.0.2 compatible)\n", SX_VERSION); return SX_TERMINATE_NOERROR; }
            const char *val = (i + 1 < argc) ? argv[i + 1] : NULL;
            if (strcmp(arg, "--mode") == 0 || strcmp(arg, "--threads") == 0) {
                if (!val) ARG_ERR();
                ++i;
                if (arg[2] == 'm') {
                    if (strcasecmp(val, "compat") == 0) s->mode = SX_MODE_COMPAT;
                    else if (strcasecmp(val, "fast") == 0) s->mode = SX_MODE_FAST;
                    else ARG_ERR();
                } else PARSE_UINT(s->threads);
                continue;
            }
            ARG_UNK();
        }
        char code = arg[1], sub_code = arg[2];
        if (toupper((unsigned char)code) == 'H' && sub_code == '\0') { sx_print_usage(); return SX_TERMINATE_NOERROR; }
        if (toupper((unsigned char)code) == 'X' && toupper((unsigned char)sub_code) == 'V') { printf("simphy-x %s (SimPhy 1.0.2 compatible)\n", SX_VERSION); return SX_TERMINATE_NOERROR; }
        ++i;
        if (i >= argc) { print_xstring_error(argv, i, "|<– missing value for this parameter\n"); return SX_SETTINGS_ERROR; }
        const char *val = argv[i];
        if (toupper((unsigned char)code) == 'I' && sub_code == '\0') continue;   /* already applied */
        long rc = apply_option(s, arg, val, argv, i);
        if (rc != SX_NO_ERROR) return rc;
    }

    /* ---- error detection (port) ---- */
    if (s->locus_tree_str != NULL && s->species_tree_str != NULL) {
        fprintf(stderr, "\n\tERROR!!! Using -S and -L (fixed species tree and locus tree) at the same time is not allowed. Please, use only a fixed locus tree (-L) or species tree (-S). In the second option, the species tree will be converted directly into a locus tree without duplication and losses. The first option allow you to add duplications (losses does not affect to the gene tree simulation.\n");
        return SX_SETTINGS_ERROR;
    }
    if ((s->species_tree_str == NULL && s->n_istrees == 0) && (s->locus_tree_str == NULL && s->n_iltrees == 0)) {
        if ((sx_is_sampling_set(s->bds_length) || sx_is_sampling_set(s->bds_leaves)) == 0) {
            fprintf(stderr, "\n\tERROR!!! You are not using neither a fixed locus tree (-L) nor a correctly specified species tree (fixed or simulated). To correctly specify a simulated species tree, you need at least one stopping rule, either bds_length or bds_leaves (-St, -Sl)\n");
            return SX_SETTINGS_ERROR;
        }
    }
    if (((s->species_tree_str != NULL || s->n_istrees > 0) && (s->locus_tree_str == NULL && s->n_iltrees == 0)) || ((s->locus_tree_str != NULL || s->n_iltrees > 0) && (s->species_tree_str == NULL && s->n_istrees == 0))) {
        if (sx_is_sampling_set(s->sb_rate) || sx_is_sampling_set(s->sd_rate) || sx_is_sampling_set(s->bds_length) || sx_is_sampling_set(s->bds_leaves)) {
            fprintf(stderr, "\n\tWARNING!!! Using a fixed tree and species tree simulation parameters has no sense, and these simulation parameters will be ignored\n");
            sx_parse_sampling("f0", &s->sb_rate, s->table); sx_parse_sampling("f0", &s->sd_rate, s->table);
            sx_parse_sampling("f0", &s->bds_length, s->table); sx_parse_sampling("f0", &s->bds_leaves, s->table);
        }
        if (s->locus_tree_str != NULL || s->n_iltrees > 0) {
            if (sx_is_variable(s->b_rate) || sx_is_variable(s->d_rate) || sx_is_variable(s->t_rate) || sx_is_variable(s->gc_rate) || s->min_lleaves > 1 || s->min_lsleaves > 1) {
                fprintf(stderr, "\n\tWARNING!!! Using locus tree birth-death parameters with a fixed locus tree has no sense, and these parameters (birth, death and minimum number of leaves) will be ignored\n");
                sx_parse_sampling("f0", &s->b_rate, s->table); sx_parse_sampling("f0", &s->d_rate, s->table);
                sx_parse_sampling("f0", &s->t_rate, s->table); sx_parse_sampling("f0", &s->gc_rate, s->table);
                s->min_lleaves = 0; s->min_lsleaves = 0;
            }
            if (sx_is_sampling_set(s->alpha_s)) {
                fprintf(stderr, "\n\tWARNING!!! Using the lineage specific substitution rate heterogeneity alpha with a fixed locus tree has no sense, and it will be ignored\n");
                sx_parse_sampling("f0", &s->alpha_s, s->table);
            }
            fprintf(stderr, "\n\tWARNING!!! Using a fixed user–specified locus tree is an advanced option, and it should not be used by general users. Some locus–tree–specific parameters may not be checked and could induce biologically senseless simulation scenarios. Examples: Different generation time for paralogs in the same species tree branch.\n");
        }
    } else if (s->species_tree_str == NULL && s->n_istrees == 0) {
        if (sx_is_sampling_set(s->bds_length) && s->min_lleaves > 2)
            fprintf(stderr, "\n\tWARNING!!! Using the minimum number of locus trees (-Ll) and a time-guided simulation of the species tree (-St) could produce incompatible values (minimum locus tree leaves > species tree leaves). If this happens, the number of minimum locus tree leaves will be changed to the number of species tree leaves.\n");
    }
    if (sx_is_sampling_set(s->gc_rate) && sx_is_sampling_set(s->b_rate) == 0) {
        fprintf(stderr, "\n\tWARNING!!! Using locus tree gene conversion simulation with no paralogs (without births) has no sense, and this parameter will be ignored\n");
        sx_parse_sampling("f0", &s->gc_rate, s->table);
    }
    fflush(stdout); fflush(stderr);
    return SX_NO_ERROR;
}

void sx_print_global_settings(FILE *file, const sx_settings *s)
{
    char buffer[256];
    const sx_stable tab = s->table;
    fprintf(file, "\nGlobal settings:\n--------------------\n\nTrees:");
    if (s->locus_tree_str != NULL || s->n_iltrees > 0) {
        fprintf(file, "\n\t-Species trees: Not taken into consideration\n");
        if (s->ltree_iname != NULL) fprintf(file, "\t-Locus trees: %d fixed trees found in %s", s->n_iltrees, s->ltree_iname);
        else fprintf(file, "\t-Locus tree: Fixed\n\t\t--> %s", s->locus_tree_str);
    } else {
        if (s->species_tree_str != NULL || s->n_istrees > 0) {
            if (s->stree_iname != NULL) fprintf(file, "\n\t-Species trees: %d fixed trees found in %s", s->n_istrees, s->stree_iname);
            else fprintf(file, "\n\t-Species tree: Fixed\n\t\t--> %s", s->species_tree_str);
        } else {
            sx_print_sampling(&s->sb_rate, buffer, tab);
            fprintf(file, "\n\t-Species trees: %d Birth-death simulations\n\t\t-Speciation rate: %s", s->ns_trees, buffer);
            sx_print_sampling(&s->sd_rate, buffer, tab);
            fprintf(file, "\n\t\t-Extinction rate: %s", buffer);
            if (sx_get(s->outgroup) > 0 || s->outgroup.distribution_code != 0) {
                sx_print_sampling(&s->outgroup, buffer, tab);
                fprintf(file, "\n\t\t-Outgroup addition:\n\t\t\t-Ingroup divergence/divergence to the ingroup ratio: %s", buffer);
            } else fprintf(file, "\n\t\t-Outgroup addition: No addition");
            if (sx_get(s->bds_leaves) > 0 || s->bds_leaves.distribution_code != 0) sx_print_sampling(&s->bds_leaves, buffer, tab);
            else sprintf(buffer, "%s", "Not used");
            fprintf(file, "\n\t\t-Stopping rules:\n\t\t\t-Number of leaves: %s", buffer);
            if (sx_get(s->bds_length) > 0 || s->bds_length.distribution_code != 0) sx_print_sampling(&s->bds_length, buffer, tab);
            else sprintf(buffer, "%s", "Not used");
            fprintf(file, "\n\t\t\t-Generations: %s", buffer);
        }
        if (!(sx_is_variable(s->lb_rate) || sx_is_variable(s->ld_rate) || sx_is_variable(s->lt_rate) || sx_is_variable(s->lgc_rate))) {
            sx_print_sampling(&s->nl_trees, buffer, tab);
            fprintf(file, "\n\t-Locus trees: %s directly obtained from each species tree (no birth-death process)\n", buffer);
        } else {
            sx_print_sampling(&s->nl_trees, buffer, tab);
            fprintf(file, "\n\t-Locus trees: %s locus simulations per species tree \n\t\t", buffer);
            sx_print_sampling(&s->b_rate, buffer, tab); fprintf(file, "-Duplication: Genome-wide parameter: %s", buffer);
            sx_print_sampling(&s->lb_rate, buffer, tab); fprintf(file, ", Rate: %s", buffer);
            sx_print_sampling(&s->d_rate, buffer, tab); fprintf(file, "\n\t\t-Loss: Genome-wide parameter: %s", buffer);
            sx_print_sampling(&s->ld_rate, buffer, tab); fprintf(file, ", Rate: %s", buffer);
            sx_print_sampling(&s->t_rate, buffer, tab); fprintf(file, "\n\t\t-Transfer: Genome-wide parameter: %s", buffer);
            sx_print_sampling(&s->lt_rate, buffer, tab); fprintf(file, ", Rate: %s", buffer);
            sx_print_sampling(&s->gc_rate, buffer, tab); fprintf(file, "\n\t\t-Gene conversion: Genome-wide parameter: %s", buffer);
            sx_print_sampling(&s->lgc_rate, buffer, tab); fprintf(file, ", Rate: %s", buffer);
            fprintf(file, "\n\t\t-Transfer (HGT and GC) sampling strategy: %s", s->t_kind == 0 ? "Random" : "Inversely proportional to the distance");
            fprintf(file, "\n\t\t-Minimum number of leaves: %u\n\t\t-Minimum number of leaves from different species: %u\n", (unsigned)s->min_lleaves, (unsigned)s->min_lsleaves);
        }
    }
    fprintf(file, "\n\t-Gene trees: %d multilocus coalescent simulations\n", s->ng_trees);
    sx_print_sampling(&s->Ne, buffer, tab); fprintf(file, "\nParameters:\n\t-Haploid efective population size: %s", buffer);
    sx_print_sampling(&s->gen_time, buffer, tab); fprintf(file, "\n\t-Generation time: %s", buffer);
    sx_print_sampling(&s->mu, buffer, tab); fprintf(file, "\n\t-Global substitution rate: %s", buffer);
    sx_print_sampling(&s->alpha_s, buffer, tab);
    fprintf(file, "\n\t-Substitution rate heterogeneities\n\t\t-Lineage (species) specific rate heterogeneity gamma shape: %s", (sx_get(s->alpha_s) == 0 && sx_is_variable(s->alpha_s) == 0) ? "No heterogeneity" : buffer);
    sx_print_sampling(&s->alpha_l, buffer, tab);
    fprintf(file, "\n\t\t-Gene family (locus tree) specific rate heterogeneity gamma shape: %s", (sx_get(s->alpha_l) == 0 && sx_is_variable(s->alpha_l) == 0) ? "No heterogeneity" : buffer);
    sx_print_sampling(&s->salpha_g, buffer, tab);
    fprintf(file, "\n\t\t-Gene tree branch specific rate heterogeneity gamma shape: Genome-wide parameter (hyperhyperparameter) %s", sx_is_sampling_set(s->salpha_g) ? buffer : "No heterogeneity");
    sx_print_sampling(&s->lalpha_g, buffer, tab);
    fprintf(file, ", locus-tree related parameter %s", (sx_get(s->lalpha_g) != 0 || sx_is_variable(s->lalpha_g) != 0) ? buffer : "No heterogeneity");
    sx_print_sampling(&s->alpha_g, buffer, tab);
    fprintf(file, ", Gamma parameter(gene-tree related) %s", (sx_get(s->alpha_g) != 0 || sx_is_variable(s->alpha_g) != 0) ? buffer : "No heterogeneity");
    sx_print_sampling(&s->ind_per_sp, buffer, tab); fprintf(file, "\n\t-Individuals per species: %s", buffer);
    fprintf(file, "\n\nMisc parameters:\n\t-Rooting method epsilon: %f\n\t-Seed: %lu", (double)s->epsilon_brent, s->u_seed);
    fprintf(file, "\n\nI/O options:\n\t-Output files prefix: %s\n\t-Verbosity: %d\n\t-Stats file: %s\n\t-Mapping: %s\n\t-Database: %s\n\t-Parameterization: %s \n\t-Command-line arguments: %s\n\t-Bounded locus subtrees: %s\n\t-Output trees with internal node labels: %s\n", s->out_name, s->verbosity, s->stats == 0 ? "OFF" : "ON", s->map == 0 ? "OFF" : "ON", s->db == 0 ? "OFF" : "ON", s->params == 0 ? "OFF" : "ON", s->command == 0 ? "OFF" : "ON", s->out_daughters == 0 ? "OFF" : "ON", s->out_labels == 0 ? "OFF" : "ON");
    fflush(file);
}

long sx_check_sampled_sloop(const sx_settings *s, int min_lleaves)
{
    int is_error = 0;
    if (sx_is_sampling_set(s->bds_leaves) && (sx_get(s->bds_leaves) <= 0)) { fprintf(stderr, "\n\tImproper value sampling the parameter -Sl, Number of species tree leaves. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_is_sampling_set(s->bds_length) && (sx_get(s->bds_length) <= 1)) { fprintf(stderr, "\n\tImproper value sampling the parameter -St, Species tree length. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_is_sampling_set(s->sb_rate) && (sx_get(s->sb_rate) <= 0)) { fprintf(stderr, "\n\tImproper value sampling the parameter -Sb, Species tree speciation rate. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_is_sampling_set(s->sd_rate) && (sx_get(s->sd_rate) > sx_get(s->sb_rate))) { fprintf(stderr, "\n\tImproper value sampling the parameter -Sd, Species tree extinction rate. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_is_sampling_set(s->outgroup) && (sx_get(s->outgroup) <= 0)) { fprintf(stderr, "\n\tImproper value sampling the parameter -So, Ingroup divergence/divergence to the ingroup ratio. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_get(s->ind_per_sp) < 1) { fprintf(stderr, "\n\tImproper value sampling the parameter -SI, Number of individuals per species. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_get(s->nl_trees) < 1) { fprintf(stderr, "\n\tImproper value sampling the parameter -Rl, Number of locus trees per species tree. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_get(s->Ne) < 2) { fprintf(stderr, "\n\tImproper value sampling the parameter -SP, Haploid effective population size. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_is_sampling_set(s->alpha_s) && (sx_get(s->alpha_s) < 0)) { fprintf(stderr, "\n\tImproper value sampling the parameter -Hs, Species-specific heterogeneity parameter. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_is_sampling_set(s->alpha_l) && (sx_get(s->alpha_l) < 0)) { fprintf(stderr, "\n\tImproper value sampling the parameter -Hl, Gene-family-specific heterogeneity parameter. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_is_sampling_set(s->salpha_g) && (sx_get(s->salpha_g) < 0)) { fprintf(stderr, "\n\tImproper value sampling the parameter -GP, Genome-wide gene-tree-branch-specific heterogeneity parameter. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_is_sampling_set(s->mu) && (sx_get(s->mu) <= 0)) { fprintf(stderr, "\n\tImproper value sampling the parameter -SU, Substitution rate. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_is_sampling_set(s->gen_time) && (sx_get(s->gen_time) <= 0)) { fprintf(stderr, "\n\tImproper value sampling the parameter -SG, Generation time. Please, check your sampling settings and try again\n"); is_error = 1; }
    if ((sx_get(s->sb_rate) < sx_get(s->sd_rate)) && (sx_get(s->bds_leaves) > 0)) { fprintf(stderr, "\n\tERROR!!! The BDSA algorithm of the conditioned birth death simulation process of the species tree does not allow birth rates less than the death rate\n"); is_error = 1; }
    else if ((sx_get(s->bds_leaves) > 0) && (sx_get(s->bds_leaves) < min_lleaves)) { fprintf(stderr, "\n\tERROR:The minimum number of locus tree nodes is bigger than the number of desired nodes of the species trees. Please, check the -Ll and -Sl parameters.\n"); is_error = 1; }
    if (is_error) { fflush(stderr); return SX_SETTINGS_ERROR; }
    return SX_NO_ERROR;
}

long sx_check_sampled_lloop(const sx_settings *s)
{
    int is_error = 0;
    if (sx_is_sampling_set(s->lb_rate) && (sx_get(s->lb_rate) < 0)) { fprintf(stderr, "\n\tImproper value sampling the parameter -Lb, Locus tree duplication rate. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_is_sampling_set(s->ld_rate) && (sx_get(s->ld_rate) > sx_get(s->lb_rate))) { fprintf(stderr, "\n\tImproper value sampling the parameter -Ld, Locus tree loss rate. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_is_sampling_set(s->lt_rate) && (sx_get(s->lt_rate) < 0)) { fprintf(stderr, "\n\tImproper value sampling the parameter -Lt, Locus tree transfer rate. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_is_sampling_set(s->lgc_rate) && (sx_get(s->lgc_rate) < 0)) { fprintf(stderr, "\n\tImproper value sampling the parameter -Lg, Locus tree gene conversion rate. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (sx_is_sampling_set(s->lalpha_g) && (sx_get(s->lalpha_g) < 0)) { fprintf(stderr, "\n\tImproper value sampling the parameter -Hh,  Gene-tree-branch-specific heterogeneity locus tree parameter. Please, check your sampling settings and try again\n"); is_error = 1; }
    if (is_error) { fflush(stderr); return SX_SETTINGS_ERROR; }
    return SX_NO_ERROR;
}

long sx_check_sampled_gloop(const sx_settings *s)
{
    if (sx_is_sampling_set(s->alpha_g) && (sx_get(s->alpha_g) < 0)) {
        fprintf(stderr, "\n\tImproper value sampling the parameter -HG, Gene-tree-branch-specific heterogeneity parameter. Please, check your sampling settings and try again\n");
        fflush(stderr);
        return SX_SETTINGS_ERROR;
    }
    return SX_NO_ERROR;
}

static sx_sunit *remap_unit(const sx_settings *src, sx_settings *dst, sx_sunit *p)
{
    const char *b = (const char *)src, *e = b + sizeof(*src), *q = (const char *)p;
    if (q >= b && q < e) return (sx_sunit *)((char *)dst + (q - b));
    return p;
}

sx_settings *sx_settings_snapshot(const sx_settings *src)
{
    sx_settings *d = sx_malloc(sizeof(*d));
    memcpy(d, src, sizeof(*d));
    for (int i = 0; i < 23; ++i) {
        d->table_storage[i].p = remap_unit(src, d, src->table_storage[i].p);
        sx_sunit *u = d->table_storage[i].p;
        for (int k = 0; k < 5; ++k)
            if (u->params_type[k] == SX_SU) u->params[k].p = remap_unit(src, d, u->params[k].p);
    }
    d->table.table = d->table_storage;
    return d;
}

void sx_settings_snapshot_free(sx_settings *s) { free(s); }

/* PrintSettingsSloop port (verbosity > 2, once per replicate). */
void sx_print_settings_sloop(const sx_settings *s, int n_replicate)
{
    char buffer[128];
    sx_out("\n\nReplicate %u parameters:\n------------------------\n\nTrees:", (unsigned)n_replicate);
    if (s->locus_tree_str != NULL || s->n_iltrees > 0) {
        sx_out("\n\t-Species tree: Not taken into consideration\n");
        if (s->ltree_iname != NULL) sx_out("\t-Locus trees: %d fixed trees found in %s", s->n_iltrees, s->ltree_iname);
        else sx_out("\t-Locus tree: Fixed\n\t\t--> %s", s->locus_tree_str);
        if (sx_get(s->gen_time) != 1) sx_out("\n\t\t-Time units: time\n\t\t-Generation time %e\n", sx_get(s->gen_time));
        else sx_out("\n\t\t-Time units: generations\n");
    } else {
        if (s->species_tree_str != NULL) {
            sx_out("\n\t-Species tree: Fixed\n\t\t--> %s", s->species_tree_str);
            if (sx_get(s->gen_time) != 1) sx_out("\n\t\t-Time units: time\n\t\t-Generation time %e\n", sx_get(s->gen_time));
            else sx_out("\n\t\t-Time units: generations\n");
        } else {
            if (sx_get(s->gen_time) != 1) {
                sx_out("\n\t-Species tree: Birth-death simulation\n\t\t-Speciation rate (per time unit): %e", sx_get(s->sb_rate));
                sx_out("\n\t\t-Extinction rate (per time unit): %e", sx_get(s->sd_rate));
                if (sx_get(s->outgroup) > 0 || s->outgroup.distribution_code != 0) sx_out("\n\t\t-Outgroup addition:\n\t\t\t-Ingroup divergence/divergence to the ingroup ratio: %e", sx_get(s->outgroup));
                else sx_out("\n\t\t-Outgroup addition: No addition");
                if (sx_get(s->bds_leaves) > 0 || s->bds_leaves.distribution_code != 0) sx_out("\n\t\t-Stopping rules:\n\t\t\t-Number of leaves: %lf", sx_get(s->bds_leaves));
                else sx_out("\n\t\t-Stopping rules:\n\t\t\t-Number of leaves: %s", "Not used");
                if (sx_get(s->bds_length) > 0 || s->bds_length.distribution_code != 0) sx_out("\n\t\t\t-Time: %e", sx_get(s->bds_length));
                else sx_out("\n\t\t\t-Time: %s", "Not used");
                sx_out("\n\t\t-Time units: time\n\t\t-Generation time %e\n", sx_get(s->gen_time));
                sx_out("\t\tNOTE: The intensity of the b-d process corresponds to the time units instead of the number of generations\n");
            } else {
                sx_out("\n\t-Species tree: Birth-death simulation\n\t\t-Speciation rate: %e", sx_get(s->sb_rate));
                sx_out("\n\t\t-Extinction rate: %e", sx_get(s->sd_rate));
                if (sx_get(s->outgroup) > 0 || s->outgroup.distribution_code != 0) sx_out("\n\t\t-Outgroup addition:\n\t\t\t-Ingroup divergence/divergence to the ingroup ratio: %e", sx_get(s->outgroup));
                else sx_out("\n\t\t-Outgroup addition: No addition");
                if (sx_get(s->bds_leaves) > 0 || s->bds_leaves.distribution_code != 0) sx_out("\n\t\t-Stopping rules:\n\t\t\t-Number of leaves: %lf", sx_get(s->bds_leaves));
                else sx_out("\n\t\t-Stopping rules:\n\t\t\t-Number of leaves: %s", "Not used");
                if (sx_get(s->bds_length) > 0 || s->bds_length.distribution_code != 0) sx_out("\n\t\t\t-Generations: %e", sx_get(s->bds_length));
                else sx_out("\n\t\t\t-Generations: %s", "Not used");
                sx_out("\n\t\t-Time units: generations\n");
            }
        }
        if (!(sx_is_variable(s->lb_rate) || sx_is_variable(s->ld_rate) || sx_is_variable(s->lt_rate) || sx_is_variable(s->lgc_rate)))
            sx_out("\n\t-%u Locus trees: directly obtained from the species tree (no birth-death process)\n", (unsigned)sx_get(s->nl_trees));
        else {
            sx_out("\n\t-%u Locus trees: birth-death simulation\n\t\t-Genome-wide duplication parameter: %e", (unsigned)sx_get(s->nl_trees), sx_get(s->b_rate));
            sx_out("\n\t\t-Genome-wide loss parameter: %e", sx_get(s->d_rate));
            sx_out("\n\t\t-Genome-wide transfer parameter: %e", sx_get(s->t_rate));
            sx_out("\n\t\t-Genome-wide gene conversion parameter: %e", sx_get(s->gc_rate));
            sx_out("\n\t\t-Minimum number of leaves: %u\n\t\t-Minimum number of leaves from different species: %u\n", (unsigned)s->min_lleaves, (unsigned)s->min_lsleaves);
        }
    }
    sx_out("\n\t-Gene trees: %d multilocus coalescent simulations\n", s->ng_trees);
    sx_out("\nParameters:\n\t-Haploid efective population size: %e", sx_get(s->Ne));
    sx_out("\n\t-Global substitution rate: %e", sx_get(s->mu));
    snprintf(buffer, sizeof buffer, "%e", sx_get(s->alpha_s));
    sx_out("\n\t-Substitution rate heterogeneities\n\t\t-Lineage (species) specific rate heterogeneity gamma shape: %s", sx_get(s->alpha_s) == 0 ? "No heterogeneity" : buffer);
    snprintf(buffer, sizeof buffer, "%e", sx_get(s->alpha_l));
    sx_out("\n\t\t-Gene family (locus tree) specific rate heterogeneity gamma shape: %s", sx_get(s->alpha_l) == 0 ? "No heterogeneity" : buffer);
    snprintf(buffer, sizeof buffer, "%e", sx_get(s->salpha_g));
    sx_out("\n\t\t-Gene tree branch specific rate heterogeneity gamma shape: Genome-wide parameter %s", sx_get(s->salpha_g) == 0 ? "No heterogeneity" : buffer);
    sx_out("\n\t-Individuals per species: %e \n\nSimulation:\n-----------\n", sx_get(s->ind_per_sp));
}

void sx_print_settings_lloop(const sx_settings *s, int curr_ltree, double gamma_l)
{
    if (sx_get(s->lb_rate) != 0 || sx_get(s->ld_rate) != 0 || sx_get(s->lt_rate) != 0 || sx_get(s->lgc_rate) != 0)
        sx_out("\n\tLocus tree %d parameters:\n\t-----------------------\n\t-Duplication rate (per generation): %e\n\t-Loss rate (per generation): %e\n\t-Transfer rate (per generation): %e\n\t-Gene conversion rate (per generation): %e\n\t-Gene family rate multiplier %e\n\t-Gene-tree-branch-specific heterogeneity locus-tree related parameter (hyperparameter) %e\n\n", curr_ltree, sx_get(s->lb_rate), sx_get(s->ld_rate), sx_get(s->lt_rate), sx_get(s->lgc_rate), gamma_l, sx_get(s->lalpha_g));
    else
        sx_out("\n\tLocus tree %d parameters:\n\t-----------------------\n\t-Fixed\n\t-Gene family rate multiplier %e\n\t-Gene-tree-branch-specific heterogeneity locus-tree related parameter (hyperparameter) %e\n\n", curr_ltree, gamma_l, sx_get(s->lalpha_g));
}

void sx_print_settings_gloop(const sx_settings *s, int curr_gtree)
{
    sx_out("\n\t\tGene tree %d parameters:\n\t\t--------------------\n\t\t-Gene tree branch specific heterogeneity parameter %e\n", curr_gtree, sx_get(s->alpha_g));
}
