/* SPDX-License-Identifier: GPL-3.0-only */
/* Compares SimPhy's double-precision lineage-count probabilities (as logged by
 * a compat run with SIMPHY_X_PROBLOG) against exact values (2048-bit MPFR). */
#include "../src/sx_coal.h"
#include "../src/sx_unif.h"
#include <mpfr.h>

static double ref_prob(int i, int o, double t, double N)
{
    if (i == o) return exp(-(double)o * (o - 1) / (2.0 * N) * t);
    int prec = 2048;
    mpfr_t sum, term, c, e;
    mpfr_inits2(prec, sum, term, c, e, (mpfr_ptr)0);
    mpfr_set_zero(sum, 1); mpfr_set_si(c, 1, MPFR_RNDN);
    for (int j = 0; j < o; ++j) { mpfr_mul_si(c, c, (long)(o + j) * (long)(i - j), MPFR_RNDN); mpfr_div_si(c, c, i + j, MPFR_RNDN); }
    mpfr_fac_ui(e, (unsigned long)o, MPFR_RNDN); mpfr_div(c, c, e, MPFR_RNDN);
    for (int k = o; k <= i; ++k) {
        if (k > o) { mpfr_mul_si(c, c, (long)(2 * k - 1) * (long)(o + k - 2), MPFR_RNDN); mpfr_mul_si(c, c, i - k + 1, MPFR_RNDN); mpfr_div_si(c, c, (long)(2 * k - 3) * (long)(k - o), MPFR_RNDN); mpfr_div_si(c, c, i + k - 1, MPFR_RNDN); }
        mpfr_set_si(e, -(long)k * (long)(k - 1), MPFR_RNDN); mpfr_mul_d(e, e, t, MPFR_RNDN); mpfr_div_d(e, e, 2.0 * N, MPFR_RNDN); mpfr_exp(e, e, MPFR_RNDN);
        mpfr_mul(term, e, c, MPFR_RNDN);
        if ((k - o) % 2 == 0) mpfr_add(sum, sum, term, MPFR_RNDN); else mpfr_sub(sum, sum, term, MPFR_RNDN);
    }
    double r = mpfr_get_d(sum, MPFR_RNDN);
    mpfr_clears(sum, term, c, e, (mpfr_ptr)0);
    return r;
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s problog\n", argv[0]); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("open"); return 2; }
    char where[8]; int i, o, Ne; double t, p;
    long n = 0, zero_but_true = 0, bad = 0, n_S = 0, bad_S = 0, zero_S = 0;
    double worst = 0; long worst_line = 0;
    double sum_ref_S = 0, sum_sp_S = 0;
    while (fscanf(f, "%7s %d %d %lf %d %lf", where, &i, &o, &t, &Ne, &p) == 6) {
        ++n;
        double ref = ref_prob(i, o, t, (double)Ne);
        double rel = ref > 0 ? fabs(p - ref) / ref : (p == 0 ? 0 : INFINITY);
        int isS = where[0] == 'S';
        if (isS) { ++n_S; sum_ref_S += ref; sum_sp_S += p; }
        if (p == 0 && ref > 1e-30) { ++zero_but_true; if (isS) ++zero_S; }
        else if (ref > 1e-30 && rel > 1e-6) { ++bad; if (isS) ++bad_S; }
        if (ref > 1e-30 && rel > worst) { worst = rel; worst_line = n; }
        if (n <= 12 || (ref > 1e-30 && (p == 0 || rel > 1e-3) && bad + zero_but_true <= 12))
            printf("%s i=%d o=%d t/N=%.3g  simphy=%.6e  exact=%.6e  rel=%.2e\n", where, i, o, t / Ne, p, ref, rel);
    }
    printf("\n%ld evaluations (C = CalcProbs, S = count sampling)\n", n);
    printf("SimPhy returned 0 while exact > 1e-30: %ld (in sampling: %ld)\n", zero_but_true, zero_S);
    printf("relative error > 1e-6 (exact > 1e-30): %ld (in sampling: %ld)\n", bad, bad_S);
    printf("worst relative error %.3g at line %ld\n", worst, worst_line);
    printf("count-sampling evaluations: %ld; sum(simphy)/sum(exact) = %.6f\n", n_S, sum_ref_S > 0 ? sum_sp_S / sum_ref_S : 0.0);
    return 0;
}
