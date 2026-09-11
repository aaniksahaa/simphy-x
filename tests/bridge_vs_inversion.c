/* Distributional comparison of the two bounded-coalescent time samplers:
 * SimPhy's sequential CDF inversion (double attempt + 512-bit MPFR redo, the
 * exact compat logic) versus the uniformisation bridge sampler.  For each
 * parameter set both samplers draw n paths; the KS statistic between the two
 * samples of every waiting time and of the total is reported. */
#include "../src/sx_unif.h"
#include "../src/sx_coal.h"

static int cmpd(const void *a, const void *b) { double x = *(const double *)a, y = *(const double *)b; return (x < y) ? -1 : (x > y); }

/* two-sample KS statistic for sorted arrays */
static double ks_stat(double *a, int na, double *b, int nb)
{
    qsort(a, (size_t)na, sizeof(double), cmpd); qsort(b, (size_t)nb, sizeof(double), cmpd);
    int i = 0, j = 0; double d = 0;
    while (i < na && j < nb) {
        if (a[i] <= b[j]) ++i; else ++j;
        double diff = fabs((double)i / na - (double)j / nb);
        if (diff > d) d = diff;
    }
    return d;
}

/* SimPhy inline logic for one bounded node */
static void simphy_waits(int i, int o, double T, int N, float eps, gsl_rng *r, double *w)
{
    int n_coals = i - o, avail = i, j;
    double cur = T; /* current_ngen relative: bound = cur - 0 */
    for (j = 0; j < n_coals; ++j) {
        w[j] = sx_sample_coal_time_mlc(avail, o, cur, N, eps, gsl_rng_uniform_pos(r), 0);
        cur -= w[j];
        if (cur - eps <= 0) break;
        --avail;
    }
    if (w[n_coals - 1] <= 0 || cur - eps <= 0 || j < n_coals) {
        cur = T; avail = i;
        for (j = 0; j < n_coals; ++j) {
            w[j] = sx_mpfr_sample_coal_time_mlc(avail, o, cur, N, eps, gsl_rng_uniform_pos(r), 0, 512);
            cur -= w[j]; --avail;
        }
    }
}

int main(int argc, char **argv)
{
    int n = argc > 1 ? atoi(argv[1]) : 20000;
    struct { int i, o; double tn; int N; } P[] = { { 40, 20, 1e-3, 290000000 }, { 63, 28, 4e-3, 290000000 }, { 30, 1, 2.0, 290000000 }, { 12, 3, 0.05, 290000000 }, { 8, 1, 0.5, 10000 }, { 80, 3, 8.0, 290000000 } };
    gsl_rng *r1 = sx_rng_fast_new(11), *r2 = sx_rng_fast_new(22);
    sx_unif_ws ws; sx_unif_ws_init(&ws);
    int worst_fail = 0;
    for (size_t c = 0; c < sizeof P / sizeof *P; ++c) {
        int i = P[c].i, o = P[c].o, N = P[c].N; double T = P[c].tn * N;
        int nc = i - o;
        double *A = sx_malloc(sizeof(double) * (size_t)n * (size_t)(nc + 1)), *B = sx_malloc(sizeof(double) * (size_t)n * (size_t)(nc + 1));
        double *w = sx_malloc(sizeof(double) * (size_t)nc);
        for (int s = 0; s < n; ++s) {
            simphy_waits(i, o, T, N, 1e-6f, r1, w);
            double tot = 0; for (int j = 0; j < nc; ++j) { A[(size_t)j * n + s] = w[j]; tot += w[j]; } A[(size_t)nc * n + s] = tot;
            sx_sample_bounded_waits(&ws, i, o, T, (double)N, 1e-6f, r2, w);
            tot = 0; for (int j = 0; j < nc; ++j) { B[(size_t)j * n + s] = w[j]; tot += w[j]; } B[(size_t)nc * n + s] = tot;
        }
        /* critical KS at alpha=0.001 for two samples of size n: 1.95*sqrt(2/n) */
        double crit = 1.95 * sqrt(2.0 / n);
        double dmax = 0; int jmax = -1;
        for (int j = 0; j <= nc; ++j) {
            double d = ks_stat(A + (size_t)j * n, n, B + (size_t)j * n, n);
            if (d > dmax) { dmax = d; jmax = j; }
        }
        double mA = 0, mB = 0; for (int s = 0; s < n; ++s) { mA += A[(size_t)nc * n + s]; mB += B[(size_t)nc * n + s]; }
        printf("i=%d o=%d T/N=%g: max KS over %d waits+total = %.4f (crit %.4f) at %s%d; mean total: simphy %.4g bridge %.4g %s\n", i, o, P[c].tn, nc + 1, dmax, crit, jmax == nc ? "total" : "wait ", jmax == nc ? 0 : jmax + 1, mA / n / N, mB / n / N, dmax > crit ? "  <-- FAIL" : "ok");
        if (dmax > crit) worst_fail = 1;
        free(A); free(B); free(w);
    }
    printf("bridges %lu fallbacks %lu\n", ws.n_bridge, ws.n_bridge_fallback);
    return worst_fail;
}
