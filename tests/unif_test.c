/* Unit test: uniformisation vs 512-bit MPFR Tavaré probabilities, and the
 * bridge sampler's first-waiting-time distribution vs the MPFR CDF. */
#include "../src/sx_unif.h"
#include "../src/sx_coal.h"
#include <mpfr.h>

extern double sx_tavare_prob(sx_unif_ws *ws, int i, int o, double t, double N);

/* reference P_{i,o}(t) at 512 bits (Tavaré) */
static double ref_prob(int i, int o, double t, double N)
{
    int prec = 2048;
    mpfr_t sum, term, c, e;
    mpfr_inits2(prec, sum, term, c, e, (mpfr_ptr)0);
    mpfr_set_zero(sum, 1);
    mpfr_set_si(c, 1, MPFR_RNDN);
    for (int j = 0; j < o; ++j) { mpfr_mul_si(c, c, (long)(o + j) * (long)(i - j), MPFR_RNDN); mpfr_div_si(c, c, i + j, MPFR_RNDN); }
    mpfr_fac_ui(e, (unsigned long)o, MPFR_RNDN); mpfr_div(c, c, e, MPFR_RNDN);
    for (int k = o; k <= i; ++k) {
        if (k > o) {
            mpfr_mul_si(c, c, (long)(2 * k - 1) * (long)(o + k - 2), MPFR_RNDN); mpfr_mul_si(c, c, i - k + 1, MPFR_RNDN);
            mpfr_div_si(c, c, (long)(2 * k - 3) * (long)(k - o), MPFR_RNDN); mpfr_div_si(c, c, i + k - 1, MPFR_RNDN);
        }
        mpfr_set_si(e, -(long)k * (long)(k - 1), MPFR_RNDN); mpfr_mul_d(e, e, t, MPFR_RNDN); mpfr_div_d(e, e, 2.0 * N, MPFR_RNDN); mpfr_exp(e, e, MPFR_RNDN);
        mpfr_mul(term, e, c, MPFR_RNDN);
        if ((k - o) % 2 == 0) mpfr_add(sum, sum, term, MPFR_RNDN); else mpfr_sub(sum, sum, term, MPFR_RNDN);
    }
    double r = mpfr_get_d(sum, MPFR_RNDN);
    mpfr_clears(sum, term, c, e, (mpfr_ptr)0);
    return r;
}

static int cmp_test(const void *x, const void *y) { double a = *(const double *)x, b = *(const double *)y; return (a < b) ? -1 : (a > b); }

int main(void)
{
    sx_unif_ws ws; sx_unif_ws_init(&ws);
    int fails = 0, n = 0;
    double worst = 0;
    int is[] = { 2, 5, 10, 25, 40, 63, 120 };
    int os[] = { 1, 2, 5, 10, 20, 30, 60 };
    double tn[] = { 1e-6, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 5.0 };
    double N = 290000000.0;
    for (size_t a = 0; a < sizeof is / sizeof *is; ++a) {
        for (size_t b = 0; b < sizeof os / sizeof *os; ++b) {
            for (size_t c = 0; c < sizeof tn / sizeof *tn; ++c) {
                int i = is[a], o = os[b];
                if (o > i) continue;
                double t = tn[c] * N;
                double ref = ref_prob(i, o, t, N);
                if (ref < 1e-300) continue;
                double p = sx_prob_coal(&ws, i, o, t, N);
                double rel = fabs(p - ref) / ref;
                ++n;
                if (rel > worst) worst = rel;
                if (rel > 1e-9 && ref > 1e-290) { ++fails; printf("FAIL i=%d o=%d t/N=%g ref=%.15g unif=%.15g rel=%.2e\n", i, o, tn[c], ref, p, rel); }
            }
        }
    }
    printf("column/prob checks: %d, failures: %d, worst rel err %.2e\n", n, fails, worst);

    /* forward transport check: v = e_j → out[k] should equal P_{j,k}(t) */
    int fmax = 40; double t = 0.02 * N;
    double *v = sx_calloc((size_t)fmax + 1, sizeof(double)), *out = sx_calloc((size_t)fmax + 1, sizeof(double));
    int ff = 0; double wf = 0;
    for (int j = 1; j <= fmax; j += 13) {
        for (int k = 0; k <= fmax; ++k) v[k] = 0;
        v[j] = 1;
        sx_unif_forward(&ws, fmax, t, N, v, out);
        for (int k = 1; k <= j; ++k) {
            double ref = ref_prob(j, k, t, N);
            if (ref < 1e-200) continue;
            double rel = fabs(out[k] - ref) / ref;
            if (rel > wf) wf = rel;
            if (rel > 1e-9) { ++ff; printf("FAIL forward j=%d k=%d ref=%.15g got=%.15g\n", j, k, ref, out[k]); }
        }
    }
    printf("forward checks worst rel err %.2e, failures %d\n", wf, ff);

    /* bridge sampler: empirical CDF of the first waiting time vs MPFR-derived CDF
     * F(w) = sum_{k<i} P_{i,k}(w) P_{k,o}(T-w) / P_{i,o}(T) */
    gsl_rng *r = sx_rng_fast_new(12345);
    int i = 30, o = 12; double T = 1e-3 * N;
    int nsamp = 200000;
    double *waits = sx_calloc((size_t)(i - o), sizeof(double));
    double qs[] = { 0.1, 0.25, 0.5, 0.75, 0.9 };
    double *first = sx_calloc((size_t)nsamp, sizeof(double));
    for (int s = 0; s < nsamp; ++s) { sx_sample_bounded_waits(&ws, i, o, T, N, 1e-6f, r, waits); first[s] = waits[0]; }
    /* sort */
    qsort(first, (size_t)nsamp, sizeof(double), cmp_test);
    double PT = ref_prob(i, o, T, N);
    int fb = 0;
    for (size_t q = 0; q < sizeof qs / sizeof *qs; ++q) {
        double w = first[(int)(qs[q] * nsamp)];
        double F = 0;
        for (int k = o; k < i; ++k) F += ref_prob(i, k, w, N) * ref_prob(k, o, T - w, N);
        F /= PT;
        double se = sqrt(qs[q] * (1 - qs[q]) / nsamp);
        printf("bridge quantile %.2f: empirical w=%.6g, theoretical CDF at w = %.5f (tolerance ~%.4f)%s\n", qs[q], w, F, 4 * se, fabs(F - qs[q]) > 4 * se ? "  FAIL" : "");
        if (fabs(F - qs[q]) > 4 * se) ++fb;
    }
    printf("bridge stats: %lu bridges, %lu fallbacks, %lu mpfr escalations\n", ws.n_bridge, ws.n_bridge_fallback, ws.n_mpfr);
    printf("%s\n", (fails || ff || fb) ? "UNIF TEST FAILED" : "UNIF TEST PASSED");
    return (fails || ff || fb) ? 1 : 0;
}
