/* One-sample tests of both bounded-coalescent samplers against the EXACT
 * marginal CDF of the j-th coalescence time τ_j given i→o lineages in T:
 *   P(τ_j <= x) = sum_{k <= i-j} P_{i,k}(x) P_{k,o}(T-x) / P_{i,o}(T),
 * with all P computed by the Tavaré series at 2048 bits. */
#include "../src/sx_unif.h"
#include "../src/sx_coal.h"
#include <mpfr.h>

static void ref_prob_mp(mpfr_t out, int i, int o, double t, double N, int prec)
{
    mpfr_t sum, term, c, e;
    mpfr_inits2(prec, sum, term, c, e, (mpfr_ptr)0);
    mpfr_set_zero(sum, 1); mpfr_set_si(c, 1, MPFR_RNDN);
    if (i < o) { mpfr_set(out, sum, MPFR_RNDN); mpfr_clears(sum, term, c, e, (mpfr_ptr)0); return; }
    for (int j = 0; j < o; ++j) { mpfr_mul_si(c, c, (long)(o + j) * (long)(i - j), MPFR_RNDN); mpfr_div_si(c, c, i + j, MPFR_RNDN); }
    mpfr_fac_ui(e, (unsigned long)o, MPFR_RNDN); mpfr_div(c, c, e, MPFR_RNDN);
    for (int k = o; k <= i; ++k) {
        if (k > o) { mpfr_mul_si(c, c, (long)(2 * k - 1) * (long)(o + k - 2), MPFR_RNDN); mpfr_mul_si(c, c, i - k + 1, MPFR_RNDN); mpfr_div_si(c, c, (long)(2 * k - 3) * (long)(k - o), MPFR_RNDN); mpfr_div_si(c, c, i + k - 1, MPFR_RNDN); }
        mpfr_set_si(e, -(long)k * (long)(k - 1), MPFR_RNDN); mpfr_mul_d(e, e, t, MPFR_RNDN); mpfr_div_d(e, e, 2.0 * N, MPFR_RNDN); mpfr_exp(e, e, MPFR_RNDN);
        mpfr_mul(term, e, c, MPFR_RNDN);
        if ((k - o) % 2 == 0) mpfr_add(sum, sum, term, MPFR_RNDN); else mpfr_sub(sum, sum, term, MPFR_RNDN);
    }
    mpfr_set(out, sum, MPFR_RNDN);
    mpfr_clears(sum, term, c, e, (mpfr_ptr)0);
}

/* exact P(τ_j <= x) */
static double exact_cdf(int i, int o, double T, double N, int j, double x)
{
    int prec = 2048;
    mpfr_t num, den, a, b;
    mpfr_inits2(prec, num, den, a, b, (mpfr_ptr)0);
    mpfr_set_zero(num, 1);
    for (int k = o; k <= i - j; ++k) {
        ref_prob_mp(a, i, k, x, N, prec);
        ref_prob_mp(b, k, o, T - x, N, prec);
        mpfr_mul(a, a, b, MPFR_RNDN);
        mpfr_add(num, num, a, MPFR_RNDN);
    }
    ref_prob_mp(den, i, o, T, N, prec);
    mpfr_div(num, num, den, MPFR_RNDN);
    double r = mpfr_get_d(num, MPFR_RNDN);
    mpfr_clears(num, den, a, b, (mpfr_ptr)0);
    return r;
}

static int cmpd(const void *a, const void *b) { double x = *(const double *)a, y = *(const double *)b; return (x < y) ? -1 : (x > y); }

static void simphy_waits(int i, int o, double T, int N, float eps, gsl_rng *r, double *w, int *used_mpfr)
{
    int n_coals = i - o, avail = i, j;
    double cur = T;
    *used_mpfr = 0;
    for (j = 0; j < n_coals; ++j) {
        w[j] = sx_sample_coal_time_mlc(avail, o, cur, N, eps, gsl_rng_uniform_pos(r), 0);
        cur -= w[j];
        if (cur - eps <= 0) break;
        --avail;
    }
    if (w[n_coals - 1] <= 0 || cur - eps <= 0 || j < n_coals) {
        *used_mpfr = 1;
        cur = T; avail = i;
        for (j = 0; j < n_coals; ++j) {
            w[j] = sx_mpfr_sample_coal_time_mlc(avail, o, cur, N, eps, gsl_rng_uniform_pos(r), 0, 512);
            cur -= w[j]; --avail;
        }
    }
}

/* one-sample KS statistic of sorted sample vs exact CDF evaluated at the sample points (subsampled for speed) */
static double ks_one(double *x, int n, int i, int o, double T, double N, int j, int npts)
{
    qsort(x, (size_t)n, sizeof(double), cmpd);
    double d = 0;
    for (int q = 0; q < npts; ++q) {
        int idx = (int)((q + 0.5) / npts * n);
        double F = exact_cdf(i, o, T, N, j, x[idx]);
        double emp = (idx + 0.5) / n;
        if (fabs(F - emp) > d) d = fabs(F - emp);
    }
    return d;
}

int main(int argc, char **argv)
{
    int i = argc > 1 ? atoi(argv[1]) : 40, o = argc > 2 ? atoi(argv[2]) : 20;
    double tn = argc > 3 ? atof(argv[3]) : 1e-3;
    int nb = argc > 4 ? atoi(argv[4]) : 20000, ns = argc > 5 ? atoi(argv[5]) : 1500;
    double N = 290000000.0, T = tn * N;
    int nc = i - o;
    gsl_rng *r1 = sx_rng_fast_new(5), *r2 = sx_rng_fast_new(6);
    sx_unif_ws ws; sx_unif_ws_init(&ws);
    double *w = sx_malloc(sizeof(double) * (size_t)nc);
    double *B = sx_malloc(sizeof(double) * (size_t)nb * (size_t)nc), *A = sx_malloc(sizeof(double) * (size_t)ns * (size_t)nc);
    for (int s = 0; s < nb; ++s) {
        sx_sample_bounded_waits(&ws, i, o, T, N, 1e-6f, r2, w);
        double acc = 0; for (int j = 0; j < nc; ++j) { acc += w[j]; B[(size_t)j * nb + s] = acc; }
    }
    int nm = 0;
    int *um_flag = sx_calloc((size_t)ns, sizeof(int));
    for (int s = 0; s < ns; ++s) {
        int um; simphy_waits(i, o, T, (int)N, 1e-6f, r1, w, &um); nm += um; um_flag[s] = um;
        double acc = 0; for (int j = 0; j < nc; ++j) { acc += w[j]; A[(size_t)j * ns + s] = acc; }
    }
    printf("i=%d o=%d T/N=%g: SimPhy used the MPFR redo in %d of %d paths\n", i, o, tn, nm, ns);
    /* split SimPhy's paths by the path that produced them */
    double *AM = sx_malloc(sizeof(double) * (size_t)ns * (size_t)nc), *AD = sx_malloc(sizeof(double) * (size_t)ns * (size_t)nc);
    int nM = 0, nD = 0;
    for (int s = 0; s < ns; ++s) {
        if (um_flag[s]) { for (int j = 0; j < nc; ++j) AM[(size_t)j * ns + nM] = A[(size_t)j * ns + s]; ++nM; }
        else { for (int j = 0; j < nc; ++j) AD[(size_t)j * ns + nD] = A[(size_t)j * ns + s]; ++nD; }
    }
    {
        int j = nc / 2;
        double *tmpM = sx_malloc(sizeof(double) * (size_t)nM), *tmpD = sx_malloc(sizeof(double) * (size_t)(nD > 0 ? nD : 1));
        for (int s = 0; s < nM; ++s) tmpM[s] = AM[(size_t)(j - 1) * ns + s];
        for (int s = 0; s < nD; ++s) tmpD[s] = AD[(size_t)(j - 1) * ns + s];
        double dM = ks_one(tmpM, nM, i, o, T, N, j, 40);
        printf("  tau_%d, SimPhy MPFR-path subset (n=%d): KS=%.4f (crit %.4f)%s\n", j, nM, dM, 1.95 / sqrt(nM), dM > 1.95 / sqrt(nM) ? " FAIL" : " ok");
        if (nD > 5) { double dD = ks_one(tmpD, nD, i, o, T, N, j, 20); printf("  tau_%d, SimPhy double-path subset (n=%d): KS=%.4f (crit %.4f)%s\n", j, nD, dD, 1.95 / sqrt(nD), dD > 1.95 / sqrt(nD) ? " FAIL" : " ok"); }
        free(tmpM); free(tmpD);
    }
    int js[] = { 1, nc / 4, nc / 2, 14 < nc ? 14 : nc - 1, nc };
    for (size_t q = 0; q < sizeof js / sizeof *js; ++q) {
        int j = js[q]; if (j < 1 || j > nc) continue;
        double db = ks_one(B + (size_t)(j - 1) * nb, nb, i, o, T, N, j, 40);
        double da = ks_one(A + (size_t)(j - 1) * ns, ns, i, o, T, N, j, 40);
        printf("  tau_%-2d: bridge KS=%.4f (crit %.4f)%s   simphy KS=%.4f (crit %.4f)%s\n", j, db, 1.95 / sqrt(nb), db > 1.95 / sqrt(nb) ? " FAIL" : " ok", da, 1.95 / sqrt(ns), da > 1.95 / sqrt(ns) ? " FAIL" : " ok");
    }
    return 0;
}
