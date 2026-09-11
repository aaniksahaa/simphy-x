/* Settings and command-line parsing: port of SimPhy's GetSettings plus the
 * simphy-x specific options (-XM mode, -XT threads). */
#ifndef SX_CLI_H
#define SX_CLI_H
#include "sx_common.h"
#include "sx_sampling.h"

typedef struct sx_settings {
    int ns_trees, ng_trees, verbosity, min_lleaves, min_lsleaves, stats, map, db, params, command, t_kind, out_time, out_daughters, out_labels;
    float epsilon_brent;
    unsigned long u_seed;
    char *species_tree_str, *locus_tree_str, *out_name, *confile_name, *stree_iname, *ltree_iname;
    char *confile_text;   /* contents of the -I file (copied to <prefix>.conf) */
    FILE *stree_ifile, *ltree_ifile;
    int n_istrees, n_iltrees;
    sx_sunit Ne, bds_leaves, ind_per_sp, nl_trees, lb_rate, ld_rate, lt_rate, lgc_rate, b_rate, d_rate, t_rate, gc_rate, sb_rate, sd_rate, bds_length, outgroup, gen_time, mu, alpha_s, alpha_l, salpha_g, lalpha_g, alpha_g;
    sx_sduple table_storage[23];
    sx_stable table;
    /* simphy-x */
    int mode;      /* SX_MODE_COMPAT / SX_MODE_FAST */
    int threads;   /* 0 = auto */
    int argc; char **argv;
} sx_settings;

sx_settings *sx_settings_new(void);
void sx_settings_free(sx_settings *s);
long sx_get_settings(sx_settings *s, int argc, char **argv);
void sx_print_usage(void);
void sx_print_global_settings(FILE *f, const sx_settings *s);
void sx_print_settings_sloop(const sx_settings *s, int n_replicate);
void sx_print_settings_lloop(const sx_settings *s, int curr_ltree, double gamma_l);
void sx_print_settings_gloop(const sx_settings *s, int curr_gtree);
long sx_check_sampled_sloop(const sx_settings *s, int min_lleaves);
long sx_check_sampled_lloop(const sx_settings *s);
long sx_check_sampled_gloop(const sx_settings *s);

#define SX_VERSION "0.1.0"

#endif

/* Deep copy of a settings struct with the internal sampling-unit pointers
 * remapped to the copy (strings are shared).  Free with sx_settings_snapshot_free. */
sx_settings *sx_settings_snapshot(const sx_settings *src);
void sx_settings_snapshot_free(sx_settings *s);
