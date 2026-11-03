/* SPDX-License-Identifier: GPL-3.0-only */
/* Orchestration.
 *
 * Compat mode: SimPhy's loop is executed sequentially on the main thread
 * ("skeleton": everything that consumes random numbers).  For every locus a
 * job is produced holding the locus tree, the (uncollapsed) gene trees, the
 * MLC traces, the deferred MPFR chains and the gamma multipliers.  Worker
 * threads evaluate the chains, replay the traces, collapse, apply the
 * multipliers, compute the DB/stats measures and write the per-locus files.
 * Ordered outputs (DB gene rows, stats lines) are committed by the main
 * thread in locus order.  Output is byte-identical for any thread count. */
#include "sx_run.h"
#include "sx_rng.h"
#include "sx_stree.h"
#include "sx_ltree.h"
#include "sx_gtree.h"
#include "sx_write.h"
#include "sx_db.h"
#include "sx_tpool.h"
#include "sx_unif.h"
#include "sx_newick.h"
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <libgen.h>
#include "sx_port.h"
#include <pthread.h>
#include <time.h>

static double now_s(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }
static double T_ltree = 0, T_msc = 0, T_mlc = 0, T_wait = 0;

typedef struct sx_locus_job {
    int rep, locus;
    int n_ldigits, n_gdigits;
    sx_ltree *lt; int lt_owned;
    int ng;
    sx_gtree **gt;
    sx_mlc_trace **trace;      /* NULL for complete (MSC) gene trees */
    sx_chain_list *chains;     /* per gene tree (only when trace) */
    double **gammas;           /* per gene tree or NULL */
    int *n_extralin;
    double *alpha_g;
    int collapse, out_labels, map, db, stats;
    const sx_names *names;
    sx_names *own_names;      /* fixed locus trees: names table owned by the job */
    int verbosity, out_time, mlc;
    sx_buf log; sx_buf *gpre, *gpost;   /* SimPhy's stdout text for this locus, flushed in locus order */
    unsigned long long lid;    /* DB id of the locus tree row (SimPhy's t_n_ltree) */
    /* results */
    long double *height_cu, *length_bl, *height_bl;
    int done;
} sx_locus_job;

typedef struct {
    sx_settings *s;
    sx_tpool *pool;
    pthread_mutex_t mu; pthread_cond_t cv;
    sx_locus_job **pending; int np, cap_p, head;
    int max_inflight;
    /* per replicate committer state */
    sqlite3 *database;
    FILE *stat_outfile;
    unsigned long long t_n_ltree;
    sx_ltree_stats *st_by_locus; int cap_st;
    int curr_stree;
} sx_ctx;

static long mkdir_cd(const char *name, int verbosity, const char *what)
{
    if (mkdir(name, S_IRWXU | S_IRWXG | S_IRWXO) != 0) {
        if (errno == EEXIST) { if (verbosity > 4) printf("\n\t WARNING, The output directory %s already exists\n", what); }
        else { fprintf(stderr, "\n\t ERROR, The output folder %s, is not accesible. Errno %d\n", name, errno); return SX_IO_ERROR; }
    }
    if (chdir(name) != 0) { fprintf(stderr, "\n\t ERROR, The output folder %s, is not accesible\n", name); return SX_IO_ERROR; }
    return SX_NO_ERROR;
}

static FILE *open_w(const char *name, const char *what)
{
    FILE *f = fopen(name, "w");
    if (!f) { perror(what); sx_check(SX_IO_ERROR, NULL); }
    return f;
}

/* ---------------- worker ---------------- */
#include <stdatomic.h>
typedef struct { sx_locus_job *job; sx_ctx *ctx; atomic_int pending_chains; } job_arg;
typedef struct { job_arg *ja; int g, idx; } chain_arg;

static void job_finalize(void *vp);

static void chain_task(void *vp)
{
    chain_arg *ca = vp;
    sx_locus_job *j = ca->ja->job;
    sx_chain_run_one(&j->chains[ca->g], j->trace[ca->g], ca->idx);
    if (atomic_fetch_sub(&ca->ja->pending_chains, 1) == 1) {
        /* last chain of this locus: finalise */
        if (ca->ja->ctx->pool) sx_tpool_submit(ca->ja->ctx->pool, job_finalize, ca->ja);
        else job_finalize(ca->ja);
    }
    free(ca);
}

/* replay traces (all chains done), collapse, multipliers, measures, write per-locus files */
static void job_finalize(void *vp)
{
    job_arg *ja = vp;
    sx_locus_job *j = ja->job;
    sx_ctx *ctx = ja->ctx;
    sx_lineage_pool pool; sx_pool_init(&pool);
    sx_buf gbuf; sx_buf_init(&gbuf);
    sx_buf mbuf; sx_buf_init(&mbuf);
    char fname[64];
    for (int g = 0; g < j->ng; ++g) {
        sx_gtree *gt = j->gt[g];
        if (j->trace && j->trace[g]) {
            if (j->chains[g].n > 0) sx_check(sx_sim_mlc_replay(j->lt, gt, &pool, j->trace[g]), ": replaying a gene tree");
            if (j->collapse) sx_gtree_collapse(gt, 1);
        }
        if (j->gammas && j->gammas[g]) sx_gtree_apply_multipliers(gt, j->gammas[g]);
        sx_write_gtree(&gbuf, gt, j->names, j->out_labels);
        /* SimPhy prints after WriteGTreeFile */
        switch (j->verbosity) {
        case 0: case 1: break;
        case 2: case 3: sx_buf_puts(&j->gpost[g], "\tDone\n"); break;
        default: { sx_buf_puts(&j->gpost[g], "\n\t\tDone: "); sx_buf tb; sx_buf_init(&tb); sx_write_gtree(&tb, gt, j->names, j->out_labels); sx_buf_append_nonl(&j->gpost[g], &tb); sx_buf_free(&tb); sx_buf_puts(&j->gpost[g], "\n"); break; }
        }
        if (j->map > 0 && j->map != 2) {
            snprintf(fname, sizeof fname, "%.*dl%.*dg.maplg", j->n_ldigits, j->locus, j->n_gdigits, g + 1);
            sx_buf_clear(&mbuf);
            sx_check(sx_write_mapping_lg(&mbuf, gt, j->names), ": obtaining LG mapping");
            sx_check(sx_buf_write_file(&mbuf, fname), ": obtaining LG mapping");
        }
        if (j->db > 0 || j->stats > 0) {
            j->height_cu[g] = sx_gtree_cu_height(gt, gt->root, j->lt->Ne);
            j->length_bl[g] = sx_gtree_bl_length(gt, gt->root);
        }
        if (j->stats > 0) j->height_bl[g] = sx_gtree_bl_height(gt, gt->root);
    }
    snprintf(fname, sizeof fname, "g_trees%.*d.trees", j->n_ldigits, j->locus);
    sx_check(sx_buf_write_file(&gbuf, fname), "Error opening g_tree file: ");
    sx_buf_free(&gbuf); sx_buf_free(&mbuf); sx_pool_free(&pool);
    pthread_mutex_lock(&ctx->mu);
    j->done = 1;
    pthread_cond_broadcast(&ctx->cv);
    pthread_mutex_unlock(&ctx->mu);
    free(ja);
}

/* entry: fan out the MPFR chains as independent tasks (load balance), then finalise */
static void worker_run(void *vp)
{
    job_arg *ja = vp;
    sx_locus_job *j = ja->job;
    int nchains = 0;
    for (int g = 0; g < j->ng; ++g) if (j->trace && j->trace[g]) nchains += j->chains[g].n;
    if (nchains == 0) { job_finalize(ja); return; }
    atomic_store(&ja->pending_chains, nchains);
    for (int g = 0; g < j->ng; ++g) {
        if (!(j->trace && j->trace[g])) continue;
        for (int c = 0; c < j->chains[g].n; ++c) {
            chain_arg *ca = sx_malloc(sizeof(*ca)); ca->ja = ja; ca->g = g; ca->idx = c;
            if (ja->ctx->pool) sx_tpool_submit(ja->ctx->pool, chain_task, ca);
            else chain_task(ca);
        }
    }
}

static void job_free(sx_locus_job *j)
{
    for (int g = 0; g < j->ng; ++g) {
        if (j->gt) sx_gtree_free(j->gt[g]);
        if (j->trace && j->trace[g]) { sx_trace_free(j->trace[g]); free(j->trace[g]); }
        if (j->chains) sx_chains_free(&j->chains[g]);
        if (j->gammas) free(j->gammas[g]);
    }
    free(j->gt); free(j->trace); free(j->chains); free(j->gammas); free(j->n_extralin); free(j->alpha_g);
    sx_buf_free(&j->log);
    for (int g = 0; g < j->ng; ++g) { sx_buf_free(&j->gpre[g]); sx_buf_free(&j->gpost[g]); }
    free(j->gpre); free(j->gpost);
    free(j->height_cu); free(j->length_bl); free(j->height_bl);
    if (j->lt_owned) sx_ltree_free(j->lt);
    if (j->own_names) sx_names_free(j->own_names);
    free(j);
}

/* commit ordered outputs of a finished job */
static void job_commit(sx_ctx *ctx, sx_locus_job *j)
{
    sx_settings *s = ctx->s;
    if (j->log.len) fwrite(j->log.data, 1, j->log.len, stdout);
    for (int g = 0; g < j->ng; ++g) {
        if (j->gpre[g].len) fwrite(j->gpre[g].data, 1, j->gpre[g].len, stdout);
        if (j->gpost[g].len) fwrite(j->gpost[g].data, 1, j->gpost[g].len, stdout);
    }
    sx_ltree_stats *st = &ctx->st_by_locus[j->locus];
    long double sum_cu = 0, sum_bl = 0; unsigned long long sum_extra = 0;
    for (int g = 0; g < j->ng; ++g) {
        if (s->stats > 0) { sum_cu += j->height_cu[g]; sum_bl += j->height_bl[g]; sum_extra += (unsigned long long)j->n_extralin[g]; }
        if (s->db > 0)
            sx_check(sx_db_write_gtree(&ctx->database, g + 1, j->lid, j->locus, j->rep, j->alpha_g[g], j->lt->n_gleaves, j->n_extralin[g], (double)j->height_cu[g], (double)j->length_bl[g]), ": writting gene tree table");
    }
    if (s->stats == 1)
        fprintf(ctx->stat_outfile, "%.*d;%d;%d;%d;%d;%d;%Le;%Le;%Le\n", j->n_ldigits, j->locus, st->losses, st->dups, st->transfr, st->leaves, st->gleaves, sum_cu / s->ng_trees, sum_bl / s->ng_trees, ((long double)sum_extra) / s->ng_trees);
}

static void commit_ready(sx_ctx *ctx, int block_until_below, int flush_all)
{
    for (;;) {
        pthread_mutex_lock(&ctx->mu);
        while (ctx->head < ctx->np && ctx->pending[ctx->head]->done) {
            sx_locus_job *j = ctx->pending[ctx->head++];
            pthread_mutex_unlock(&ctx->mu);
            job_commit(ctx, j);
            job_free(j);
            pthread_mutex_lock(&ctx->mu);
        }
        int inflight = ctx->np - ctx->head;
        int need_wait = flush_all ? inflight > 0 : inflight >= block_until_below;
        pthread_mutex_unlock(&ctx->mu);
        if (!need_wait) break;
        /* help: run a queued job on this thread, else wait for a completion */
        if (ctx->pool && sx_tpool_help(ctx->pool)) continue;
        pthread_mutex_lock(&ctx->mu);
        if (ctx->head < ctx->np && !ctx->pending[ctx->head]->done) pthread_cond_wait(&ctx->cv, &ctx->mu);
        pthread_mutex_unlock(&ctx->mu);
    }
    if (flush_all) { ctx->np = 0; ctx->head = 0; }
}

static void job_submit(sx_ctx *ctx, sx_locus_job *j)
{
    pthread_mutex_lock(&ctx->mu);
    if (ctx->np == ctx->cap_p) { ctx->cap_p = ctx->cap_p ? ctx->cap_p * 2 : 256; ctx->pending = sx_realloc(ctx->pending, sizeof(*ctx->pending) * (size_t)ctx->cap_p); }
    ctx->pending[ctx->np++] = j;
    pthread_mutex_unlock(&ctx->mu);
    job_arg *ja = sx_malloc(sizeof(*ja)); ja->job = j; ja->ctx = ctx; atomic_store(&ja->pending_chains, 0);
    if (ctx->pool) sx_tpool_submit(ctx->pool, worker_run, ja);
    else worker_run(ja);
}

/* ---------------- main loop ---------------- */
long sx_run(sx_settings *s)
{
    if (s->mode == SX_MODE_FAST) return sx_run_fast(s);
    sx_verbosity = s->verbosity;
    sx_mpfr_rescue = s->rescue;
    gsl_rng *r = sx_rng_compat_new(s->u_seed);
    sx_names *names = NULL;
    sx_buf buf; sx_buf_init(&buf);
    char *out_name = s->out_name;
    int verbosity = s->verbosity;
    int min_lleaves = s->min_lleaves;

    sx_ctx ctx; memset(&ctx, 0, sizeof ctx);
    ctx.s = s;
    pthread_mutex_init(&ctx.mu, NULL); pthread_cond_init(&ctx.cv, NULL);
    int nthreads = s->threads > 0 ? s->threads : sx_hardware_threads();
    ctx.pool = nthreads > 1 ? sx_tpool_new(nthreads - 1) : NULL;
    /* the probability rows stay on the main thread: workers are saturated by the MPFR chains */
    sx_calc_probs_set_pool(NULL);
    ctx.max_inflight = 2 * nthreads + 4;

    if (out_name == NULL) { out_name = s->out_name = sx_malloc(sizeof("SimPhy_outfiles")); strcpy(out_name, "SimPhy_outfiles"); }
    if (verbosity) sx_print_global_settings(stdout, s);

    int n_sdigits = (int)sx_count_intdigits((long)s->ns_trees, 0);
    umask(0);
    if (mkdir_cd(out_name, verbosity, out_name) != SX_NO_ERROR) return SX_IO_ERROR;
    char *out_copy = sx_malloc(strlen(out_name) + 1); strcpy(out_copy, out_name);
    char *base_out_name = basename(out_copy);
    char *fname = sx_malloc(strlen(base_out_name) + 64);

    if (s->db > 0) { sprintf(fname, "%s.db", base_out_name); sx_check(sx_db_init(&ctx.database, fname), NULL); }
    if (s->params > 0) {
        sprintf(fname, "%s.params", base_out_name);
        FILE *pf = open_w(fname, "Error opening .params file:");
        sx_print_global_settings(pf, s);
        fclose(pf);
    }
    if (s->command > 0) {
        sprintf(fname, "%s.command", base_out_name);
        FILE *cf = open_w(fname, "Error opening .command file:");
        for (int i = 0; i < s->argc; ++i) fprintf(cf, "%s ", s->argv[i]);
        if (s->confile_text) {
            sprintf(fname, "%s.conf", base_out_name);
            FILE *co = open_w(fname, "Error opening output .conf file:");
            fputs(s->confile_text, co);
            fclose(co);
        }
        fclose(cf);
    }
    if (verbosity > 0 && verbosity < 3) printf("\n\nSimulation:\n-----------\n");

    sx_sunit *sloop[17] = { &s->bds_leaves, &s->bds_length, &s->sb_rate, &s->sd_rate, &s->outgroup, &s->ind_per_sp, &s->nl_trees, &s->b_rate, &s->d_rate, &s->t_rate, &s->gc_rate, &s->Ne, &s->alpha_s, &s->alpha_l, &s->salpha_g, &s->mu, &s->gen_time };
    sx_sunit *lloop[5] = { &s->lb_rate, &s->ld_rate, &s->lt_rate, &s->lgc_rate, &s->lalpha_g };
    sx_sunit *gloop[1] = { &s->alpha_g };

    int *node_ptrs = NULL;
    char *nexus_ltree_str = NULL;
    char curr_outdir[32];
    sx_lineage_pool pool; sx_pool_init(&pool);
    int locus_simulated = (sx_is_variable(s->lb_rate) || sx_is_variable(s->ld_rate) || sx_is_variable(s->lt_rate) || sx_is_variable(s->lgc_rate));

    if (s->n_istrees > 0) sx_check(sx_nexus_init(s->stree_ifile), ": parsing Nexus file");
    for (int curr_stree = 1; curr_stree <= s->ns_trees; ++curr_stree) {
        ctx.curr_stree = curr_stree;
        sx_check(sx_sample_distr(r, 17, sloop), ": sampling distributions");
        sx_check(sx_check_sampled_sloop(s, min_lleaves), ": improper sampled values");

        int n_ldigits = (int)sx_count_intdigits((long)sx_get(s->nl_trees), 0);
        int n_gdigits = (int)sx_count_intdigits((long)s->ng_trees, 0);
        int nl_trees = (int)sx_get(s->nl_trees);
        double gen_time = sx_get(s->gen_time) != 1 ? sx_get(s->gen_time) : 1;
        if (nl_trees + 1 > ctx.cap_st) { ctx.st_by_locus = sx_realloc(ctx.st_by_locus, sizeof(sx_ltree_stats) * (size_t)(nl_trees + 1)); ctx.cap_st = nl_trees + 1; }

        sprintf(curr_outdir, "%.*d", n_sdigits, curr_stree);
        if (mkdir_cd(curr_outdir, verbosity, curr_outdir) != SX_NO_ERROR) return SX_IO_ERROR;

        FILE *s_outfile = NULL, *l_outfile = NULL, *d_outfile = NULL;
        ctx.stat_outfile = NULL;
        if (s->locus_tree_str == NULL && s->n_iltrees == 0) s_outfile = open_w("s_tree.trees", "Error opening s_trees.tree file..");
        if (s->stats == 1) {
            ctx.stat_outfile = open_w("stats.txt", "Error opening stats.txt:");
            fprintf(ctx.stat_outfile, "L_tree;N_losses;N_duplications;N_transfers;N_lt_totaltips;N_gt_presenttips;Mean_gt_height(cu);Mean_gt_height(ec);Mean_extra_lineages\n");
        }
        l_outfile = open_w("l_trees.trees", "Error opening l_trees.tree file: ");
        if (s->out_daughters) d_outfile = open_w("bounded_locus_subtrees.out", "Error opening bounded_locus_subtrees.out file: ");

        /* ---- species tree ---- */
        if (s->n_istrees > 0) sx_check(sx_nexus_next(s->stree_ifile, &s->species_tree_str), " : getting next Nexus tree");
        else if (s->n_iltrees > 0) sx_check(sx_nexus_init(s->ltree_ifile), ": Parsing Nexus tree");
        if (verbosity > 2) sx_print_settings_sloop(s, curr_stree);
        sx_stree *sp_tree = NULL;
        if (names) { sx_names_free(names); names = NULL; }
        int fixed_locus = (s->locus_tree_str != NULL || s->n_iltrees > 0);
        if (fixed_locus) {
            /* no species tree at all */
        } else if (s->n_istrees > 0 || s->species_tree_str != NULL) {
            sp_tree = s->n_istrees > 0 ? sx_parse_nexus_stree(s->species_tree_str, &names, verbosity, gen_time, (int)sx_get(s->Ne), sx_get(s->mu), (int)sx_get(s->ind_per_sp))
                                       : sx_parse_newick_stree(s->species_tree_str, &names, verbosity, gen_time, (int)sx_get(s->Ne), sx_get(s->mu), (int)sx_get(s->ind_per_sp));
            if (sp_tree->n_leaves < min_lleaves)
                sx_check(SX_SETTINGS_ERROR, "\n\t\tERROR:The minimum number of locus tree leaves is bigger than the number of leaves of the fixed species tree. Please, check the -Ll parameter and the input species tree.\n");
            sx_stree_collapse(sp_tree, 0);
        } else {
            if (verbosity > 2) printf("\nSimulating the species tree...");
            sx_check(sx_new_bds_tree(&sp_tree, (int)sx_get(s->bds_leaves), sx_get(s->bds_length), sx_get(s->sb_rate), sx_get(s->sd_rate), gen_time, (int)sx_get(s->Ne), sx_get(s->mu), (int)sx_get(s->ind_per_sp), sx_get(s->outgroup), 0, 1, r, verbosity), ": simulating the species tree");
            sx_stree_collapse(sp_tree, 0);
            if (verbosity > 3) { sx_buf_clear(&buf); sx_write_stree(&buf, sp_tree, NULL, s->out_time, s->out_labels); size_t n = buf.len; while (n > 0 && buf.data[n - 1] == '\n') --n; fwrite(buf.data, 1, n, stdout); printf("\n"); }
            if (sp_tree->n_leaves < min_lleaves) {
                if (verbosity > 2) printf("\n\t\tThe number of minimum locus tree leaves will be changed to the number of species tree leaves to avoid infinite loops.\n\t\t");
                min_lleaves = sp_tree->n_leaves;
            }
        }
        if (!fixed_locus) {
            if (sx_get(s->alpha_s) != 0) {
                if (verbosity > 2) printf("\tGenerating lineage specific substitution rate heterogeneity...");
                sx_stree_rateheter(sp_tree, sx_get(s->alpha_s), r);
                if (verbosity > 2) printf(" Done\n");
            }
            if (verbosity > 2) printf("Done\n");
            if (locus_simulated && node_ptrs == NULL) node_ptrs = sx_calloc((size_t)MAX_LEAVES, sizeof(int));
        }

        switch (verbosity) {
        case 1: printf("Replicate %u of %u: Simulating %u gene trees from %d locus trees... ", curr_stree, s->ns_trees, (unsigned)(nl_trees * s->ng_trees), nl_trees); break;
        case 2: printf("\nReplicate %u of %u: Simulating %u gene trees from %d locus trees...\n", curr_stree, s->ns_trees, (unsigned)(nl_trees * s->ng_trees), nl_trees); break;
        default: break;
        }
        fflush(stdout);

        if (s_outfile) { sx_buf_clear(&buf); sx_write_stree(&buf, sp_tree, names, s->out_time, s->out_labels); sx_buf_append_file(&buf, s_outfile); }
        if (s->db > 0 && fixed_locus) {
            sx_check(sx_db_write_stree(&ctx.database, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ": writting the species tree table");
        } else if (s->db > 0) {
            long double st_height = sx_stree_cu_height(sp_tree, sp_tree->root, sp_tree->Ne);
            long double st_length = sx_stree_cu_length(sp_tree, sp_tree->root, sp_tree->Ne);
            sx_check(sx_db_write_stree(&ctx.database, sp_tree->n_leaves, sx_get(s->sb_rate), sx_get(s->sd_rate), (int)sx_get(s->bds_leaves), sx_get(s->bds_length), (double)st_height, (double)st_length, sp_tree->gen_length[SX_CHILD(sp_tree, sp_tree->root, 0)], (int)sx_get(s->ind_per_sp), (int)sx_get(s->nl_trees), sx_get(s->b_rate), sx_get(s->d_rate), sx_get(s->t_rate), sx_get(s->gc_rate), sx_get(s->alpha_s), sx_get(s->alpha_l), sx_get(s->salpha_g), (int)sx_get(s->Ne), sx_get(s->mu), sx_get(s->gen_time)), ": writting the species tree table");
        }
        /* species index labels: post-order after the first simulated locus tree
         * (SimPhy: ReindexSTree after every SimBDLHTree).  Set once here so
         * that workers see stable labels. */
        if (locus_simulated && sp_tree) sx_stree_reindex(sp_tree, 1);

        /* ---- locus loop ---- */
        sx_ltree *shared_lt = NULL;   /* MSC: locus tree copied from the species tree once */
        double gamma_l = 0;
        int collapse = (s->map > 0 || s->out_labels == 1) ? 1 : 0;
        for (int curr_ltree = 1; curr_ltree <= nl_trees; ++curr_ltree) {
            sx_check(sx_sample_distr(r, 5, lloop), ": sampling distributions");
            sx_check(sx_check_sampled_lloop(s), ": improper sampled values");
            sx_ltree_stats st; memset(&st, 0, sizeof st);
            if (sx_get(s->alpha_l) != 0) gamma_l = gsl_ran_gamma(r, sx_get(s->alpha_l), 1 / (sx_get(s->alpha_l)));
            int simulated_locus = (sx_get(s->lb_rate) != 0 || sx_get(s->ld_rate) != 0 || sx_get(s->lt_rate) != 0 || sx_get(s->lgc_rate) != 0);

            sx_locus_job *job = sx_calloc(1, sizeof(*job));
            job->rep = curr_stree; job->locus = curr_ltree; job->n_ldigits = n_ldigits; job->n_gdigits = n_gdigits;
            job->verbosity = verbosity; job->out_time = s->out_time;
            sx_buf_init(&job->log);
            job->gpre = sx_calloc((size_t)s->ng_trees, sizeof(sx_buf)); job->gpost = sx_calloc((size_t)s->ng_trees, sizeof(sx_buf));
            for (int g = 0; g < s->ng_trees; ++g) { sx_buf_init(&job->gpre[g]); sx_buf_init(&job->gpost[g]); }
            sx_out_capture(&job->log);
            if (verbosity > 2) sx_print_settings_lloop(s, curr_ltree, gamma_l);
            job->ng = s->ng_trees; job->collapse = collapse; job->out_labels = s->out_labels; job->map = s->map; job->db = s->db; job->stats = s->stats;
            job->names = names;
            job->gt = sx_calloc((size_t)job->ng, sizeof(*job->gt));
            job->n_extralin = sx_calloc((size_t)job->ng, sizeof(int));
            job->alpha_g = sx_calloc((size_t)job->ng, sizeof(double));
            job->height_cu = sx_calloc((size_t)job->ng, sizeof(long double));
            job->length_bl = sx_calloc((size_t)job->ng, sizeof(long double));
            job->height_bl = sx_calloc((size_t)job->ng, sizeof(long double));

            sx_ltree *locus_tree;
            if (simulated_locus) {
                if (verbosity > 2) sx_out("\tSimulating locus tree %d ... ", curr_ltree);
                locus_tree = NULL;
                double t0 = now_s();
                sx_check(sx_sim_bdlh_tree(sp_tree, &locus_tree, node_ptrs, sx_get(s->lb_rate), sx_get(s->ld_rate), sx_get(s->lt_rate), sx_get(s->lgc_rate), s->t_kind, r, min_lleaves, s->min_lsleaves, verbosity, &st, 1), ": simulating a locus tree");
                T_ltree += now_s() - t0;
                sx_stree_reindex(sp_tree, 1);
                if (getenv("SIMPHY_X_COMPAT_EXACTPROBS") && (st.dups > 0 || st.transfr > 0 || st.gc > 0)) {
                    /* diagnostic: exact (uniformisation) lineage-count probabilities inside the compat pipeline */
                    static __thread sx_unif_ws *dws = NULL;
                    if (!dws) { dws = sx_malloc(sizeof(*dws)); sx_unif_ws_init(dws); }
                    sx_ltree_free_probs(locus_tree);
                    sx_calc_probs_nlineages_fast(locus_tree, locus_tree->root, locus_tree->Ne, 0, dws);
                }
                if (s->map > 1) locus_tree->n_gleaves += st.losses;
                if (s->map > 0 && s->map != 2) {
                    sprintf(fname, "%.*d.mapsl", n_ldigits, curr_ltree);
                    sx_buf_clear(&buf);
                    sx_check(sx_write_mapping_sl(&buf, sp_tree, locus_tree, names), "");
                    sx_check(sx_buf_write_file(&buf, fname), "");
                }
                job->lt = locus_tree; job->lt_owned = 1;
            } else if (fixed_locus) {
                if (verbosity > 2) sx_out("\tObtaining the preset locus tree (there is no birth-death process)... ");
                names = NULL;
                int nd = 0, nl = 0, nt = 0, ngc = 0;
                if (s->n_iltrees > 0) {
                    sx_check(sx_nexus_next(s->ltree_ifile, &nexus_ltree_str), ": parsing the next Nexus tree");
                    locus_tree = sx_parse_nexus_ltree(nexus_ltree_str, &names, verbosity, gen_time, (int)sx_get(s->Ne), sx_get(s->mu), (int)sx_get(s->ind_per_sp), &nd, &nl, &nt, &ngc);
                } else
                    locus_tree = sx_parse_newick_ltree(s->locus_tree_str, &names, verbosity, gen_time, (int)sx_get(s->Ne), sx_get(s->mu), (int)sx_get(s->ind_per_sp), &nd, &nl, &nt, &ngc);
                st.dups = nd; st.losses = nl; st.transfr = nt; st.gc = ngc;
                if (s->stats == 1) { st.leaves = locus_tree->n_leaves; st.gleaves = locus_tree->n_gleaves; }
                job->lt = locus_tree; job->lt_owned = 1; job->names = names; job->own_names = names;
            } else {
                if (shared_lt == NULL) {
                    shared_lt = sx_ltree_new(sp_tree->n_nodes, sp_tree->n_leaves, sp_tree->n_gleaves, sp_tree->max_children, sp_tree->gen_time, sp_tree->Ne, sp_tree->mu);
                    if (verbosity > 2) sx_out("\tInvariable locus tree: Copying the species tree as locus tree... ");
                    sx_ltree_copy_from_stree(shared_lt, sp_tree);
                } else if (verbosity > 2) sx_out("\tInvariable locus tree: Reusing the previous locus tree... ");
                locus_tree = shared_lt;
                if (s->stats == 1) { st.leaves = locus_tree->n_leaves; st.gleaves = locus_tree->n_gleaves; }
                job->lt = locus_tree; job->lt_owned = 0;
            }
            switch (verbosity) {
            case 0: case 1: case 2: break;
            case 3: sx_out(" Done"); break;
            default: { sx_out("\n\tDone: "); sx_buf tb; sx_buf_init(&tb); sx_write_ltree(&tb, locus_tree, names, s->out_time, s->out_labels); sx_buf_append_nonl(&job->log, &tb); sx_buf_free(&tb); break; }
            }
            if (sx_get(s->alpha_l) != 0) {
                if (verbosity > 2) sx_out("\n\tGenerating gene family (locus tree) specific substitution rate heterogeneity...");
                locus_tree->mu = (sp_tree ? sp_tree->mu : sx_get(s->mu)) * gamma_l;
                if (verbosity > 2) sx_out(" Done\n");
            }
            ctx.st_by_locus[curr_ltree] = st;

            /* locus tree outputs (ordered, sequential) */
            sx_buf_clear(&buf); sx_write_ltree(&buf, locus_tree, names, s->out_time, s->out_labels); sx_buf_append_file(&buf, l_outfile);
            if (s->out_daughters) {
                sx_buf_clear(&buf);
                if (st.dups != 0 || st.transfr != 0 || st.gc != 0) sx_write_daughters(&buf, locus_tree, names);
                else sx_buf_putc(&buf, '\n');
                sx_buf_append_file(&buf, d_outfile);
            }
            if (s->db > 0) {
                sx_check(sx_db_write_ltree(&ctx.database, curr_ltree, curr_stree, sx_get(s->lb_rate), sx_get(s->ld_rate), sx_get(s->lt_rate), sx_get(s->lgc_rate), locus_tree->n_leaves, st.dups, st.losses, st.transfr, st.gc, st.ntrials, gamma_l, sx_get(s->lalpha_g)), ": writting the locus tree rable");
                job->lid = ++ctx.t_n_ltree;
            }

            /* gene trees: skeleton */
            int mlc = (st.dups != 0 || st.transfr != 0 || st.gc != 0);
            if (mlc) { job->trace = sx_calloc((size_t)job->ng, sizeof(*job->trace)); job->chains = sx_calloc((size_t)job->ng, sizeof(*job->chains)); }
            job->mlc = mlc;
            for (int g = 0; g < job->ng; ++g) {
                sx_out_capture(&job->gpre[g]);
                sx_check(sx_sample_distr(r, 1, gloop), ": sampling distributions");
                sx_check(sx_check_sampled_gloop(s), ": improper sampled values");
                job->alpha_g[g] = sx_get(s->alpha_g);
                if (verbosity > 2) sx_print_settings_gloop(s, g + 1);
                switch (verbosity) {
                case 0: case 1: break;
                case 2: sx_out("\tSimulating gene tree %d... ", (curr_ltree - 1) * s->ng_trees + g + 1); break;
                default: sx_out("\n\t\tSimulating gene tree %d using a %s process...", g + 1, mlc ? "multilocus coalescent" : "multispecies coalescent"); break;
                }
                sx_gtree *gt = sx_gtree_new((locus_tree->n_gleaves * 2) - 1, locus_tree->max_children, locus_tree->gen_time);
                job->gt[g] = gt;
                double tg0 = now_s();
                if (!mlc) {
                    sx_check(sx_sim_msc_gtree(locus_tree, gt, &pool, r, &job->n_extralin[g], s->map > 1 ? 1 : 0, verbosity, collapse), ": simulating a gene tree");
                    T_msc += now_s() - tg0;
                } else {
                    job->trace[g] = sx_calloc(1, sizeof(sx_mlc_trace)); sx_trace_init(job->trace[g]); sx_chains_init(&job->chains[g]);
                    sx_check(sx_sim_mlc_gtree(locus_tree, gt, &pool, s->epsilon_brent, r, &job->n_extralin[g], verbosity, collapse, &job->chains[g], job->trace[g]), ": simulating a gene tree");
                    if (job->chains[g].n == 0) {
                        /* nothing deferred: finish here (collapse) — identical to inline */
                        if (collapse) sx_gtree_collapse(gt, 1);
                        sx_trace_free(job->trace[g]); free(job->trace[g]); job->trace[g] = NULL;
                    }
                    T_mlc += now_s() - tg0;
                }
                if (sx_get(s->alpha_g) != 0) {
                    if (verbosity > 2) sx_out("\n\t\tGenerating gene tree branch specific substitution rate heterogeneity...");
                    if (!job->gammas) job->gammas = sx_calloc((size_t)job->ng, sizeof(*job->gammas));
                    job->gammas[g] = sx_malloc(sizeof(double) * (size_t)gt->n_nodes);
                    sx_gtree_draw_multipliers(gt->n_nodes, sx_get(s->alpha_g), r, job->gammas[g]);
                    if (verbosity > 3) sx_out(" Done");
                    else if (verbosity > 2) sx_out(" Done\n\t");
                }
            }
            sx_out_capture(NULL);
            /* hand over to the workers; commit whatever is finished, keep the window bounded */
            double tc0 = now_s();
            job_submit(&ctx, job);
            commit_ready(&ctx, ctx.max_inflight, 0);
            T_wait += now_s() - tc0;
        }
        double tf0 = now_s();
        commit_ready(&ctx, 0, 1);
        T_wait += now_s() - tf0;
        if (chdir("..") != 0) { fprintf(stderr, "\n\t ERROR, impossible to come back to the parent directory\n"); return SX_IO_ERROR; }
        if (shared_lt) sx_ltree_free(shared_lt);
        if (sp_tree) sx_stree_free(sp_tree);
        if (fixed_locus) names = NULL;
        if (l_outfile) fclose(l_outfile);
        if (d_outfile) fclose(d_outfile);
        if (ctx.stat_outfile) fclose(ctx.stat_outfile);
        if (s_outfile) fclose(s_outfile);
        if (verbosity >= 1 && verbosity <= 4) printf("Done\n");
    }
    if (verbosity == 0) printf("Done\n");
    else if (verbosity > 2) printf("\nDone \n");
    if (sx_mpfr_rescue_count) fprintf(stderr, "\n[simphy-x] note: %lu bounded-coalescent inversions needed more than 512 bits; SimPhy 1.0.2 would have aborted here (-XR 1 rescue). Output is deterministic but no longer byte-comparable to SimPhy from the first such locus on.\n", sx_mpfr_rescue_count);
    if (getenv("SIMPHY_X_DEBUG")) fprintf(stderr, "[simphy-x compat] skeleton: locus trees %.2fs, MSC %.2fs, MLC pass A %.2fs; waiting for workers/commits %.2fs\n", T_ltree, T_msc, T_mlc, T_wait);

    if (names) sx_names_free(names);
    if (s->db > 0) sx_check(sx_db_close(&ctx.database), NULL);
    sx_calc_probs_set_pool(NULL);
    if (ctx.pool) sx_tpool_free(ctx.pool);
    pthread_mutex_destroy(&ctx.mu); pthread_cond_destroy(&ctx.cv);
    free(ctx.pending); free(ctx.st_by_locus);
    free(node_ptrs); free(fname); free(out_copy); free(nexus_ltree_str);
    sx_pool_free(&pool);
    sx_buf_free(&buf);
    gsl_rng_free(r);
    fflush(stdout);
    return SX_NO_ERROR;
}
