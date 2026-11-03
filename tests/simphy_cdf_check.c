/* SPDX-License-Identifier: GPL-3.0-only */
/* Evaluates SimPhy's conditional CDF of the next coalescence time (MPFR and double
 * versions, exactly as SimPhy computes them) against the exact CDF
 *   F(w) = sum_{k<i} P_{i,k}(w) P_{k,o}(T-w) / P_{i,o}(T)  (2048-bit MPFR). */
#include "../src/sx_unif.h"
#include "../src/sx_coal.h"
#include <mpfr.h>

static void ref_prob_mp(mpfr_t out, int i, int o, double t, double N, int prec)
{
    mpfr_t sum, term, c, e;
    mpfr_inits2(prec, sum, term, c, e, (mpfr_ptr)0);
    mpfr_set_zero(sum, 1); mpfr_set_si(c, 1, MPFR_RNDN);
    if (i >= o) {
        for (int j = 0; j < o; ++j) { mpfr_mul_si(c, c, (long)(o + j) * (long)(i - j), MPFR_RNDN); mpfr_div_si(c, c, i + j, MPFR_RNDN); }
        mpfr_fac_ui(e, (unsigned long)o, MPFR_RNDN); mpfr_div(c, c, e, MPFR_RNDN);
        for (int k = o; k <= i; ++k) {
            if (k > o) { mpfr_mul_si(c, c, (long)(2 * k - 1) * (long)(o + k - 2), MPFR_RNDN); mpfr_mul_si(c, c, i - k + 1, MPFR_RNDN); mpfr_div_si(c, c, (long)(2 * k - 3) * (long)(k - o), MPFR_RNDN); mpfr_div_si(c, c, i + k - 1, MPFR_RNDN); }
            mpfr_set_si(e, -(long)k * (long)(k - 1), MPFR_RNDN); mpfr_mul_d(e, e, t, MPFR_RNDN); mpfr_div_d(e, e, 2.0 * N, MPFR_RNDN); mpfr_exp(e, e, MPFR_RNDN);
            mpfr_mul(term, e, c, MPFR_RNDN);
            if ((k - o) % 2 == 0) mpfr_add(sum, sum, term, MPFR_RNDN); else mpfr_sub(sum, sum, term, MPFR_RNDN);
        }
    }
    mpfr_set(out, sum, MPFR_RNDN);
    mpfr_clears(sum, term, c, e, (mpfr_ptr)0);
}
static double exact_cdf(int i, int o, double T, double N, double x)
{
    int prec = 2048;
    mpfr_t num, den, a, b;
    mpfr_inits2(prec, num, den, a, b, (mpfr_ptr)0);
    mpfr_set_zero(num, 1);
    for (int k = o; k <= i - 1; ++k) { ref_prob_mp(a, i, k, x, N, prec); ref_prob_mp(b, k, o, T - x, N, prec); mpfr_mul(a, a, b, MPFR_RNDN); mpfr_add(num, num, a, MPFR_RNDN); }
    ref_prob_mp(den, i, o, T, N, prec);
    mpfr_div(num, num, den, MPFR_RNDN);
    double r = mpfr_get_d(num, MPFR_RNDN);
    mpfr_clears(num, den, a, b, (mpfr_ptr)0);
    return r;
}

int main(int argc, char **argv)
{
    int i = argc > 1 ? atoi(argv[1]) : 40, o = argc > 2 ? atoi(argv[2]) : 20;
    double tn = argc > 3 ? atof(argv[3]) : 1e-3;
    double N = 290000000.0, T = tn * N;
    printf("i=%d o=%d T/N=%g\n   w/T      exact F(w)   SimPhy MPFR F(w)  SimPhy double F(w)\n", i, o, tn);
    double fr[] = { 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999999 };
    for (size_t q = 0; q < sizeof fr / sizeof *fr; ++q) {
        double w = fr[q] * T;
        printf("  %8.6f  %12.8f  %14.8f  %14.8f\n", fr[q], exact_cdf(i, o, T, N, w), sx_debug_mpfr_cdf(i, o, T, (int)N, w, 512), sx_debug_double_cdf(i, o, T, (int)N, w));
    }
    return 0;
}
