/* Fast mode orchestration.
 *
 * Every (replicate, locus) pair is an independent job with its own random
 * stream derived from the seed, so the loci of a replicate (and of successive
 * replicates) run fully in parallel.  Workers own a private copy of the
 * settings (locus-level sampling units are sampled into it), a private clone
 * of the species tree (the locus-tree simulation annotates it), thread-local
 * scratch, and write the per-locus files directly.  Ordered outputs
 * (l_trees.trees, bounded_locus_subtrees.out, DB rows, stats) are committed
 * by the main thread in (replicate, locus) order, so the output is identical
 * for any thread count. */
#include "sx_run.h"
#include "sx_rng.h"
#include "sx_stree.h"
#include "sx_ltree.h"
#include "sx_gtree.h"
#include "sx_unif.h"
#include "sx_write.h"
#include "sx_db.h"
#include "sx_tpool.h"
#include "sx_newick.h"
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <libgen.h>
#include <pthread.h>

enum { TAG_REP = 1, TAG_LOCUS = 2 };

typedef struct fast_rep {
    int rep;
    sx_settings *snap;          /* replicate-level sampled values */
    sx_stree *sp_tree;          /* read-only for workers (NULL with a fixed locus tree) */
    char **ltree_strs;          /* -LR: one Nexus tree string per locus */
    int fixed_locus;
    sx_names *names;            /* fixed species tree names or NULL */
    int nl_trees, n_ldigits, n_gdigits;
    char dir[32];
    FILE *l_outfile, *d_outfile, *stat_outfile;   /* owned by the committer */
    sx_buf log; int log_flushed;                  /* species-level stdout text, printed before the first locus commit */
    int committed;              /* loci committed so far */
    int refs;
} fast_rep;

typedef struct fast_job {
    fast_rep *rep;
    int locus;
    /* outputs */
    sx_buf ltree_line, daughters_line;
    char db_ltree_sql[1024];
    int has_db_ltree;
    int ng;
    double *alpha_g; int *n_extralin; int n_gleaves;
    long double *height_cu, *length_bl, *height_bl;
    sx_ltree_stats st;
    double gamma_l;
    sx_buf log;                 /* SimPhy's stdout text for this locus */
    int done;
} fast_job;

typedef struct {
    sx_settings *s;
    sx_tpool *pool;
    pthread_mutex_t mu; pthread_cond_t cv;
    fast_job **pending; int np, cap_p, head;
    int max_inflight;
    sqlite3 *database;
    unsigned long long t_n_ltree;
} fast_ctx;

static long mkdir_only(const char *name)
{
    if (mkdir(name, S_IRWXU | S_IRWXG | S_IRWXO) != 0 && errno != EEXIST) {
        fprintf(stderr, "\n\t ERROR, The output folder %s, is not accesible. Errno %d\n", name, errno); return SX_IO_ERROR;
    }
    return SX_NO_ERROR;
}
static FILE *open_w(const char *name, const char *what)
{
    FILE *f = fopen(name, "w");
    if (!f) { perror(what); sx_check(SX_IO_ERROR, NULL); }
    return f;
}

/* ---- thread-local scratch ---- */
static __thread gsl_rng *tl_rng = NULL;
static __thread sx_unif_ws *tl_ws = NULL;
static __thread sx_lineage_pool *tl_pool = NULL;
static __thread int *tl_node_ptrs = NULL;
static pthread_mutex_t stat_mu = PTHREAD_MUTEX_INITIALIZER;
static unsigned long g_bridge = 0, g_bridge_fb = 0, g_pieces = 0, g_mpfr = 0;

static void tl_init(void)
{
    if (!tl_rng) tl_rng = sx_rng_fast_new(1);
    if (!tl_ws) { tl_ws = sx_malloc(sizeof(*tl_ws)); sx_unif_ws_init(tl_ws); }
    if (!tl_pool) { tl_pool = sx_malloc(sizeof(*tl_pool)); sx_pool_init(tl_pool); }
    if (!tl_node_ptrs) tl_node_ptrs = sx_calloc((size_t)MAX_LEAVES, sizeof(int));
}

typedef struct { fast_job *job; fast_ctx *ctx; } job_arg;

static void fast_worker(void *vp)
{
    job_arg *ja = vp;
    fast_job *j = ja->job;
    fast_ctx *ctx = ja->ctx;
    fast_rep *R = j->rep;
    sx_settings *s = sx_settings_snapshot(R->snap);
    tl_init();
    gsl_rng *r = tl_rng;
    sx_rng_fast_reseed(r, sx_stream_seed((uint64_t)s->u_seed, (uint64_t)R->rep, (uint64_t)j->locus, TAG_LOCUS));
    sx_stree *sp = R->sp_tree ? sx_stree_clone(R->sp_tree) : NULL;
    int verbosity = 0;
    char fname[128];
    sx_buf buf; sx_buf_init(&buf);
    sx_names *names = R->names;

    sx_sunit *lloop[5] = { &s->lb_rate, &s->ld_rate, &s->lt_rate, &s->lgc_rate, &s->lalpha_g };
    sx_sunit *gloop[1] = { &s->alpha_g };
    verbosity = s->verbosity;
    sx_buf_init(&j->log);
    sx_out_capture(&j->log);
    sx_check(sx_sample_distr(r, 5, lloop), ": sampling distributions");
    sx_check(sx_check_sampled_lloop(s), ": improper sampled values");
    sx_ltree_stats st; memset(&st, 0, sizeof st);
    double gamma_l = 0;
    if (sx_get(s->alpha_l) != 0) gamma_l = gsl_ran_gamma(r, sx_get(s->alpha_l), 1 / (sx_get(s->alpha_l)));
    if (verbosity > 2) sx_print_settings_lloop(s, j->locus, gamma_l);
    int simulated_locus = (sx_get(s->lb_rate) != 0 || sx_get(s->ld_rate) != 0 || sx_get(s->lt_rate) != 0 || sx_get(s->lgc_rate) != 0);
    int min_lleaves = s->min_lleaves;
    if (sp && sp->n_leaves < min_lleaves) min_lleaves = sp->n_leaves;
    int collapse = (s->map > 0 || s->out_labels == 1) ? 1 : 0;

    sx_ltree *lt = NULL;
    sx_names *own_names = NULL;
    if (R->fixed_locus) {
        if (verbosity > 2) sx_out("\tObtaining the preset locus tree (there is no birth-death process)... ");
        int nd = 0, nl = 0, nt = 0, ngc = 0;
        double gt = sx_get(s->gen_time) != 1 ? sx_get(s->gen_time) : 1;
        if (R->ltree_strs) lt = sx_parse_nexus_ltree(R->ltree_strs[j->locus - 1], &own_names, verbosity, gt, (int)sx_get(s->Ne), sx_get(s->mu), (int)sx_get(s->ind_per_sp), &nd, &nl, &nt, &ngc);
        else lt = sx_parse_newick_ltree(s->locus_tree_str, &own_names, verbosity, gt, (int)sx_get(s->Ne), sx_get(s->mu), (int)sx_get(s->ind_per_sp), &nd, &nl, &nt, &ngc);
        names = own_names;
        st.dups = nd; st.losses = nl; st.transfr = nt; st.gc = ngc;
        if (s->stats == 1) { st.leaves = lt->n_leaves; st.gleaves = lt->n_gleaves; }
    } else if (simulated_locus) {
        if (verbosity > 2) sx_out("\tSimulating locus tree %d ... ", j->locus);
        sx_check(sx_sim_bdlh_tree(sp, &lt, tl_node_ptrs, sx_get(s->lb_rate), sx_get(s->ld_rate), sx_get(s->lt_rate), sx_get(s->lgc_rate), s->t_kind, r, min_lleaves, s->min_lsleaves, verbosity, &st, getenv("SIMPHY_X_FAST_COMPATPROBS") != NULL), ": simulating a locus tree");
        sx_stree_reindex(sp, 1);
        /* replace the compat probabilities by uniformisation ones */
        /* diagnostic: SIMPHY_X_FAST_COMPATPROBS keeps SimPhy's double-precision probabilities */
        if ((st.dups > 0 || st.transfr > 0 || st.gc > 0) && !getenv("SIMPHY_X_FAST_COMPATPROBS")) sx_calc_probs_nlineages_fast(lt, lt->root, lt->Ne, 0, tl_ws);
        if (s->map > 1) lt->n_gleaves += st.losses;
        if (s->map > 0 && s->map != 2) {
            snprintf(fname, sizeof fname, "%s/%.*d.mapsl", R->dir, R->n_ldigits, j->locus);
            sx_buf_clear(&buf);
            sx_check(sx_write_mapping_sl(&buf, sp, lt, names), "");
            sx_check(sx_buf_write_file(&buf, fname), "");
        }
    } else {
        if (verbosity > 2) sx_out(j->locus == 1 ? "\tInvariable locus tree: Copying the species tree as locus tree... " : "\tInvariable locus tree: Reusing the previous locus tree... ");
        lt = sx_ltree_new(sp->n_nodes, sp->n_leaves, sp->n_gleaves, sp->max_children, sp->gen_time, sp->Ne, sp->mu);
        sx_ltree_copy_from_stree(lt, sp);
        if (s->stats == 1) { st.leaves = lt->n_leaves; st.gleaves = lt->n_gleaves; }
    }
    switch (verbosity) {
    case 0: case 1: case 2: break;
    case 3: sx_out(" Done"); break;
    default: { sx_out("\n\tDone: "); sx_buf tb; sx_buf_init(&tb); sx_write_ltree(&tb, lt, names, s->out_time, s->out_labels); sx_buf_append_nonl(&j->log, &tb); sx_buf_free(&tb); break; }
    }
    if (sx_get(s->alpha_l) != 0) {
        if (verbosity > 2) sx_out("\n\tGenerating gene family (locus tree) specific substitution rate heterogeneity...");
        lt->mu = (sp ? sp->mu : sx_get(s->mu)) * gamma_l;
        if (verbosity > 2) sx_out(" Done\n");
    }
    j->st = st; j->gamma_l = gamma_l; j->n_gleaves = lt->n_gleaves;

    sx_write_ltree(&j->ltree_line, lt, names, s->out_time, s->out_labels);
    if (s->out_daughters) {
        if (st.dups != 0 || st.transfr != 0 || st.gc != 0) sx_write_daughters(&j->daughters_line, lt, names);
        else sx_buf_putc(&j->daughters_line, '\n');
    }
    if (s->db > 0) {
        sx_db_ltree_sql(j->db_ltree_sql, sizeof j->db_ltree_sql, j->locus, R->rep, sx_get(s->lb_rate), sx_get(s->ld_rate), sx_get(s->lt_rate), sx_get(s->lgc_rate), lt->n_leaves, st.dups, st.losses, st.transfr, st.gc, st.ntrials, gamma_l, sx_get(s->lalpha_g));
        j->has_db_ltree = 1;
    }

    int mlc = (st.dups != 0 || st.transfr != 0 || st.gc != 0);
    sx_buf gbuf; sx_buf_init(&gbuf);
    sx_gtree *gt = sx_gtree_new((lt->n_gleaves * 2) - 1, lt->max_children, lt->gen_time);
    for (int g = 0; g < j->ng; ++g) {
        sx_check(sx_sample_distr(r, 1, gloop), ": sampling distributions");
        sx_check(sx_check_sampled_gloop(s), ": improper sampled values");
        j->alpha_g[g] = sx_get(s->alpha_g);
        if (verbosity > 2) sx_print_settings_gloop(s, g + 1);
        switch (verbosity) {
        case 0: case 1: break;
        case 2: sx_out("\tSimulating gene tree %d... ", (j->locus - 1) * s->ng_trees + g + 1); break;
        default: sx_out("\n\t\tSimulating gene tree %d using a %s process...", g + 1, mlc ? "multilocus coalescent" : "multispecies coalescent"); break;
        }
        if (!mlc) sx_check(sx_sim_msc_gtree(lt, gt, tl_pool, r, &j->n_extralin[g], s->map > 1 ? 1 : 0, 0, collapse), ": simulating a gene tree");
        else sx_check(sx_sim_mlc_gtree_fast(lt, gt, tl_pool, s->epsilon_brent, r, &j->n_extralin[g], collapse, tl_ws, getenv("SIMPHY_X_FAST_COMPATPROBS") != NULL), ": simulating a gene tree");
        if (sx_get(s->alpha_g) != 0) {
            if (verbosity > 2) sx_out("\n\t\tGenerating gene tree branch specific substitution rate heterogeneity...");
            sx_gtree_rateheter(gt, sx_get(s->alpha_g), r);
            if (verbosity > 3) sx_out(" Done");
            else if (verbosity > 2) sx_out(" Done\n\t");
        }
        sx_write_gtree(&gbuf, gt, names, s->out_labels);
        switch (verbosity) {
        case 0: case 1: break;
        case 2: case 3: sx_out("\tDone\n"); break;
        default: sx_out("\n\t\tDone: "); { sx_buf tb; sx_buf_init(&tb); sx_write_gtree(&tb, gt, names, s->out_labels); sx_buf_append_nonl(&j->log, &tb); sx_buf_free(&tb); } sx_out("\n"); break;
        }
        if (s->map > 0 && s->map != 2) {
            snprintf(fname, sizeof fname, "%s/%.*dl%.*dg.maplg", R->dir, R->n_ldigits, j->locus, R->n_gdigits, g + 1);
            sx_buf_clear(&buf);
            sx_check(sx_write_mapping_lg(&buf, gt, names), ": obtaining LG mapping");
            sx_check(sx_buf_write_file(&buf, fname), ": obtaining LG mapping");
        }
        if (s->db > 0 || s->stats > 0) { j->height_cu[g] = sx_gtree_cu_height(gt, gt->root, lt->Ne); j->length_bl[g] = sx_gtree_bl_length(gt, gt->root); }
        if (s->stats > 0) j->height_bl[g] = sx_gtree_bl_height(gt, gt->root);
    }
    snprintf(fname, sizeof fname, "%s/g_trees%.*d.trees", R->dir, R->n_ldigits, j->locus);
    sx_check(sx_buf_write_file(&gbuf, fname), "Error opening g_tree file: ");
    sx_buf_free(&gbuf); sx_buf_free(&buf);
    sx_gtree_free(gt); sx_ltree_free(lt); if (sp) sx_stree_free(sp); if (own_names) sx_names_free(own_names); sx_settings_snapshot_free(s);
    sx_out_capture(NULL);
    pthread_mutex_lock(&stat_mu);
    g_bridge += tl_ws->n_bridge; g_bridge_fb += tl_ws->n_bridge_fallback; g_pieces += tl_ws->n_forward_tavare; g_mpfr += tl_ws->n_mpfr;
    tl_ws->n_bridge = tl_ws->n_bridge_fallback = tl_ws->n_forward_tavare = tl_ws->n_mpfr = 0;
    pthread_mutex_unlock(&stat_mu);

    pthread_mutex_lock(&ctx->mu);
    j->done = 1;
    pthread_cond_broadcast(&ctx->cv);
    pthread_mutex_unlock(&ctx->mu);
    free(ja);
}

static void rep_open_files(fast_ctx *ctx, fast_rep *R)
{
    sx_settings *s = ctx->s;
    char p[96];
    snprintf(p, sizeof p, "%s/l_trees.trees", R->dir); R->l_outfile = open_w(p, "Error opening l_trees.tree file: ");
    if (s->out_daughters) { snprintf(p, sizeof p, "%s/bounded_locus_subtrees.out", R->dir); R->d_outfile = open_w(p, "Error opening bounded_locus_subtrees.out file: "); }
    if (s->stats == 1) {
        snprintf(p, sizeof p, "%s/stats.txt", R->dir); R->stat_outfile = open_w(p, "Error opening stats.txt:");
        fprintf(R->stat_outfile, "L_tree;N_losses;N_duplications;N_transfers;N_lt_totaltips;N_gt_presenttips;Mean_gt_height(cu);Mean_gt_height(ec);Mean_extra_lineages\n");
    }
}

static void rep_release(fast_rep *R)
{
    if (--R->refs > 0) return;
    if (R->l_outfile) fclose(R->l_outfile);
    if (R->d_outfile) fclose(R->d_outfile);
    if (R->stat_outfile) fclose(R->stat_outfile);
    sx_buf_free(&R->log);
    if (R->ltree_strs) { for (int i = 0; i < R->nl_trees; ++i) free(R->ltree_strs[i]); free(R->ltree_strs); }
    if (R->sp_tree) sx_stree_free(R->sp_tree);
    if (R->names) sx_names_free(R->names);
    sx_settings_snapshot_free(R->snap);
    free(R);
}

static void job_commit(fast_ctx *ctx, fast_job *j)
{
    sx_settings *s = ctx->s;
    fast_rep *R = j->rep;
    if (!R->log_flushed) { if (R->log.len) fwrite(R->log.data, 1, R->log.len, stdout); R->log_flushed = 1; }
    if (j->log.len) fwrite(j->log.data, 1, j->log.len, stdout);
    sx_buf_free(&j->log);
    sx_buf_append_file(&j->ltree_line, R->l_outfile);
    if (s->out_daughters) sx_buf_append_file(&j->daughters_line, R->d_outfile);
    unsigned long long lid = 0;
    if (j->has_db_ltree) { sx_check(sx_db_exec(&ctx->database, j->db_ltree_sql), ": writting the locus tree rable"); lid = ++ctx->t_n_ltree; }
    long double sum_cu = 0, sum_bl = 0; unsigned long long sum_extra = 0;
    for (int g = 0; g < j->ng; ++g) {
        if (s->stats > 0) { sum_cu += j->height_cu[g]; sum_bl += j->height_bl[g]; sum_extra += (unsigned long long)j->n_extralin[g]; }
        if (s->db > 0) sx_check(sx_db_write_gtree(&ctx->database, g + 1, lid, j->locus, R->rep, j->alpha_g[g], j->n_gleaves, j->n_extralin[g], (double)j->height_cu[g], (double)j->length_bl[g]), ": writting gene tree table");
    }
    if (s->stats == 1)
        fprintf(R->stat_outfile, "%.*d;%d;%d;%d;%d;%d;%Le;%Le;%Le\n", R->n_ldigits, j->locus, j->st.losses, j->st.dups, j->st.transfr, j->st.leaves, j->st.gleaves, sum_cu / s->ng_trees, sum_bl / s->ng_trees, ((long double)sum_extra) / s->ng_trees);
    R->committed++;
    if (R->committed == R->nl_trees && s->verbosity >= 1 && s->verbosity <= 4) { printf("Done\n"); fflush(stdout); }
    sx_buf_free(&j->ltree_line); sx_buf_free(&j->daughters_line);
    free(j->alpha_g); free(j->n_extralin); free(j->height_cu); free(j->length_bl); free(j->height_bl);
    rep_release(R);
    free(j);
}

static void commit_ready(fast_ctx *ctx, int block_until_below, int flush_all)
{
    for (;;) {
        pthread_mutex_lock(&ctx->mu);
        while (ctx->head < ctx->np && ctx->pending[ctx->head]->done) {
            fast_job *j = ctx->pending[ctx->head++];
            pthread_mutex_unlock(&ctx->mu);
            job_commit(ctx, j);
            pthread_mutex_lock(&ctx->mu);
        }
        if (ctx->head == ctx->np) { ctx->head = ctx->np = 0; }
        int inflight = ctx->np - ctx->head;
        int need_wait = flush_all ? inflight > 0 : inflight >= block_until_below;
        pthread_mutex_unlock(&ctx->mu);
        if (!need_wait) break;
        if (ctx->pool && sx_tpool_help(ctx->pool)) continue;
        pthread_mutex_lock(&ctx->mu);
        if (ctx->head < ctx->np && !ctx->pending[ctx->head]->done) pthread_cond_wait(&ctx->cv, &ctx->mu);
        pthread_mutex_unlock(&ctx->mu);
    }
}

static void job_submit(fast_ctx *ctx, fast_job *j)
{
    pthread_mutex_lock(&ctx->mu);
    if (ctx->np == ctx->cap_p) { ctx->cap_p = ctx->cap_p ? ctx->cap_p * 2 : 256; ctx->pending = sx_realloc(ctx->pending, sizeof(*ctx->pending) * (size_t)ctx->cap_p); }
    ctx->pending[ctx->np++] = j;
    pthread_mutex_unlock(&ctx->mu);
    job_arg *ja = sx_malloc(sizeof(*ja)); ja->job = j; ja->ctx = ctx;
    if (ctx->pool) sx_tpool_submit(ctx->pool, fast_worker, ja);
    else fast_worker(ja);
}

long sx_run_fast(sx_settings *s)
{
    sx_verbosity = s->verbosity;
    int verbosity = s->verbosity;
    fast_ctx ctx; memset(&ctx, 0, sizeof ctx);
    ctx.s = s;
    pthread_mutex_init(&ctx.mu, NULL); pthread_cond_init(&ctx.cv, NULL);
    int nthreads = s->threads > 0 ? s->threads : sx_hardware_threads();
    ctx.pool = nthreads > 1 ? sx_tpool_new(nthreads - 1) : NULL;
    ctx.max_inflight = 4 * nthreads + 8;
    char *out_name = s->out_name;
    if (out_name == NULL) { out_name = s->out_name = sx_malloc(sizeof("SimPhy_outfiles")); strcpy(out_name, "SimPhy_outfiles"); }
    if (verbosity) sx_print_global_settings(stdout, s);

    int n_sdigits = (int)sx_count_intdigits((long)s->ns_trees, 0);
    umask(0);
    if (mkdir_only(out_name) != SX_NO_ERROR) return SX_IO_ERROR;
    if (chdir(out_name) != 0) { fprintf(stderr, "\n\t ERROR, The output folder %s, is not accesible\n", out_name); return SX_IO_ERROR; }
    char *out_copy = sx_malloc(strlen(out_name) + 1); strcpy(out_copy, out_name);
    char *base_out_name = basename(out_copy);
    char *fname = sx_malloc(strlen(base_out_name) + 64);
    if (s->db > 0) { sprintf(fname, "%s.db", base_out_name); sx_check(sx_db_init(&ctx.database, fname), NULL); }
    if (s->params > 0) { sprintf(fname, "%s.params", base_out_name); FILE *pf = open_w(fname, "Error opening .params file:"); sx_print_global_settings(pf, s); fclose(pf); }
    if (s->command > 0) { sprintf(fname, "%s.command", base_out_name); FILE *cf = open_w(fname, "Error opening .command file:"); for (int i = 0; i < s->argc; ++i) fprintf(cf, "%s ", s->argv[i]);
        if (s->confile_text) {
            sprintf(fname, "%s.conf", base_out_name);
            FILE *co = open_w(fname, "Error opening output .conf file:");
            fputs(s->confile_text, co);
            fclose(co);
        } fclose(cf); }
    if (verbosity > 0 && verbosity < 3) printf("\n\nSimulation:\n-----------\n");

    sx_sunit *sloop[17] = { &s->bds_leaves, &s->bds_length, &s->sb_rate, &s->sd_rate, &s->outgroup, &s->ind_per_sp, &s->nl_trees, &s->b_rate, &s->d_rate, &s->t_rate, &s->gc_rate, &s->Ne, &s->alpha_s, &s->alpha_l, &s->salpha_g, &s->mu, &s->gen_time };
    int locus_simulated = (sx_is_variable(s->lb_rate) || sx_is_variable(s->ld_rate) || sx_is_variable(s->lt_rate) || sx_is_variable(s->lgc_rate));
    gsl_rng *r = sx_rng_fast_new(1);
    sx_buf buf; sx_buf_init(&buf);
    int min_lleaves = s->min_lleaves;

    if (s->n_istrees > 0) sx_check(sx_nexus_init(s->stree_ifile), ": parsing Nexus file");
    for (int curr_stree = 1; curr_stree <= s->ns_trees; ++curr_stree) {
        sx_rng_fast_reseed(r, sx_stream_seed((uint64_t)s->u_seed, (uint64_t)curr_stree, 0, TAG_REP));
        sx_check(sx_sample_distr(r, 17, sloop), ": sampling distributions");
        sx_check(sx_check_sampled_sloop(s, min_lleaves), ": improper sampled values");
        fast_rep *R = sx_calloc(1, sizeof(*R));
        R->rep = curr_stree;
        sx_buf_init(&R->log);
        sx_out_capture(&R->log);
        if (s->n_istrees > 0) sx_check(sx_nexus_next(s->stree_ifile, &s->species_tree_str), " : getting next Nexus tree");
        else if (s->n_iltrees > 0) sx_check(sx_nexus_init(s->ltree_ifile), ": Parsing Nexus tree");
        if (verbosity > 2) sx_print_settings_sloop(s, curr_stree);
        R->nl_trees = (int)sx_get(s->nl_trees);
        if (s->n_iltrees > 0) {
            R->ltree_strs = sx_calloc((size_t)R->nl_trees, sizeof(char *));
            for (int l = 0; l < R->nl_trees; ++l) sx_check(sx_nexus_next(s->ltree_ifile, &R->ltree_strs[l]), ": parsing the next Nexus tree");
        }
        R->n_ldigits = (int)sx_count_intdigits((long)R->nl_trees, 0);
        R->n_gdigits = (int)sx_count_intdigits((long)s->ng_trees, 0);
        double gen_time = sx_get(s->gen_time) != 1 ? sx_get(s->gen_time) : 1;
        snprintf(R->dir, sizeof R->dir, "%.*d", n_sdigits, curr_stree);
        if (mkdir_only(R->dir) != SX_NO_ERROR) return SX_IO_ERROR;

        sx_stree *sp_tree = NULL;
        sx_names *names = NULL;
        R->fixed_locus = (s->locus_tree_str != NULL || s->n_iltrees > 0);
        if (R->fixed_locus) {
            /* no species tree */
        } else if (s->n_istrees > 0 || s->species_tree_str != NULL) {
            sp_tree = s->n_istrees > 0 ? sx_parse_nexus_stree(s->species_tree_str, &names, verbosity, gen_time, (int)sx_get(s->Ne), sx_get(s->mu), (int)sx_get(s->ind_per_sp))
                                       : sx_parse_newick_stree(s->species_tree_str, &names, verbosity, gen_time, (int)sx_get(s->Ne), sx_get(s->mu), (int)sx_get(s->ind_per_sp));
            if (sp_tree->n_leaves < min_lleaves)
                sx_check(SX_SETTINGS_ERROR, "\n\t\tERROR:The minimum number of locus tree leaves is bigger than the number of leaves of the fixed species tree. Please, check the -Ll parameter and the input species tree.\n");
            sx_stree_collapse(sp_tree, 0);
        } else {
            if (verbosity > 2) sx_out("\nSimulating the species tree...");
            sx_check(sx_new_bds_tree(&sp_tree, (int)sx_get(s->bds_leaves), sx_get(s->bds_length), sx_get(s->sb_rate), sx_get(s->sd_rate), gen_time, (int)sx_get(s->Ne), sx_get(s->mu), (int)sx_get(s->ind_per_sp), sx_get(s->outgroup), 0, 1, r, verbosity), ": simulating the species tree");
            sx_stree_collapse(sp_tree, 0);
            if (verbosity > 3) { sx_buf tb; sx_buf_init(&tb); sx_write_stree(&tb, sp_tree, NULL, s->out_time, s->out_labels); sx_buf_append_nonl(&R->log, &tb); sx_buf_free(&tb); sx_out("\n"); }
            if (sp_tree->n_leaves < min_lleaves) {
                if (verbosity > 2) sx_out("\n\t\tThe number of minimum locus tree leaves will be changed to the number of species tree leaves to avoid infinite loops.\n\t\t");
                min_lleaves = sp_tree->n_leaves;
            }
        }
        if (!R->fixed_locus) {
            if (sx_get(s->alpha_s) != 0) {
                if (verbosity > 2) sx_out("\tGenerating lineage specific substitution rate heterogeneity...");
                sx_stree_rateheter(sp_tree, sx_get(s->alpha_s), r);
                if (verbosity > 2) sx_out(" Done\n");
            }
            if (verbosity > 2) sx_out("Done\n");
        }
        switch (verbosity) {
        case 1: sx_out("Replicate %u of %u: Simulating %u gene trees from %d locus trees... ", curr_stree, s->ns_trees, (unsigned)(R->nl_trees * s->ng_trees), R->nl_trees); break;
        case 2: sx_out("\nReplicate %u of %u: Simulating %u gene trees from %d locus trees...\n", curr_stree, s->ns_trees, (unsigned)(R->nl_trees * s->ng_trees), R->nl_trees); break;
        default: break;
        }
        sx_out_capture(NULL);

        if (!R->fixed_locus) {
            sprintf(fname, "%s/s_tree.trees", R->dir);
            sx_buf_clear(&buf); sx_write_stree(&buf, sp_tree, names, s->out_time, s->out_labels);
            sx_check(sx_buf_write_file(&buf, fname), "Error opening s_trees.tree file..");
        }
        if (s->db > 0 && R->fixed_locus) {
            sx_check(sx_db_write_stree(&ctx.database, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ": writting the species tree table");
        } else if (s->db > 0) {
            long double st_height = sx_stree_cu_height(sp_tree, sp_tree->root, sp_tree->Ne);
            long double st_length = sx_stree_cu_length(sp_tree, sp_tree->root, sp_tree->Ne);
            sx_check(sx_db_write_stree(&ctx.database, sp_tree->n_leaves, sx_get(s->sb_rate), sx_get(s->sd_rate), (int)sx_get(s->bds_leaves), sx_get(s->bds_length), (double)st_height, (double)st_length, sp_tree->gen_length[SX_CHILD(sp_tree, sp_tree->root, 0)], (int)sx_get(s->ind_per_sp), (int)sx_get(s->nl_trees), sx_get(s->b_rate), sx_get(s->d_rate), sx_get(s->t_rate), sx_get(s->gc_rate), sx_get(s->alpha_s), sx_get(s->alpha_l), sx_get(s->salpha_g), (int)sx_get(s->Ne), sx_get(s->mu), sx_get(s->gen_time)), ": writting the species tree table");
        }
        if (locus_simulated && sp_tree) sx_stree_reindex(sp_tree, 1);
        R->sp_tree = sp_tree; R->names = names;
        R->snap = sx_settings_snapshot(s);
        R->snap->min_lleaves = min_lleaves;
        rep_open_files(&ctx, R);
        R->refs = R->nl_trees;
        for (int curr_ltree = 1; curr_ltree <= R->nl_trees; ++curr_ltree) {
            fast_job *j = sx_calloc(1, sizeof(*j));
            j->rep = R; j->locus = curr_ltree; j->ng = s->ng_trees;
            sx_buf_init(&j->ltree_line); sx_buf_init(&j->daughters_line);
            j->alpha_g = sx_calloc((size_t)j->ng, sizeof(double)); j->n_extralin = sx_calloc((size_t)j->ng, sizeof(int));
            j->height_cu = sx_calloc((size_t)j->ng, sizeof(long double)); j->length_bl = sx_calloc((size_t)j->ng, sizeof(long double)); j->height_bl = sx_calloc((size_t)j->ng, sizeof(long double));
            job_submit(&ctx, j);
            commit_ready(&ctx, ctx.max_inflight, 0);
        }
    }
    commit_ready(&ctx, 0, 1);
    if (verbosity == 0) printf("Done\n");
    else if (verbosity > 2) printf("\nDone \n");
    if (getenv("SIMPHY_X_DEBUG")) fprintf(stderr, "[simphy-x] bridges=%lu bridge_fallbacks=%lu forward_pieces=%lu mpfr_escalations=%lu\n", g_bridge, g_bridge_fb, g_pieces, g_mpfr);
    if (s->db > 0) sx_check(sx_db_close(&ctx.database), NULL);
    if (ctx.pool) sx_tpool_free(ctx.pool);
    pthread_mutex_destroy(&ctx.mu); pthread_cond_destroy(&ctx.cv);
    free(ctx.pending); free(fname); free(out_copy);
    sx_buf_free(&buf);
    gsl_rng_free(r);
    fflush(stdout);
    return SX_NO_ERROR;
}
