/* SPDX-License-Identifier: GPL-3.0-only */
/* Fast-mode numerics for the bounded coalescent.
 *
 * Inside a bounded locus subtree the lineage count along a branch is a pure
 * death chain with rates lam_k = k(k-1)/(2N).  All quantities are computed by
 * uniformisation (positive series, no cancellation, double precision):
 *   - forward transport of a lineage-count distribution through a branch
 *     (replaces CalcProbsNLineagesLTree's O(n^3) alternating series),
 *   - the column P_{k,o}(t) for all k (lineage-count sampling),
 *   - exact sampling of endpoint-conditioned coalescence times with the
 *     Hobolth–Stone (2009) bridge construction (replaces CDF inversion in MPFR).
 * For very large Lambda*t (many expected events) the well-conditioned Tavaré
 * series is used, with adaptive-precision MPFR as the last resort. */
#ifndef SX_UNIF_H
#define SX_UNIF_H
#include "sx_common.h"
#include "sx_rng.h"

/* Thread-local scratch (grown on demand). */
typedef struct sx_unif_ws {
    double *w; int w_cap;          /* Poisson weights */
    double *b; size_t b_cap;       /* backward vectors */
    double *v0, *v1; int v_cap;    /* running vectors */
    double *tau; int tau_cap;      /* event times */
    unsigned long n_bridge, n_bridge_fallback, n_forward_tavare, n_mpfr;   /* statistics */
} sx_unif_ws;

void sx_unif_ws_init(sx_unif_ws *ws);
void sx_unif_ws_free(sx_unif_ws *ws);

/* Threshold on Lambda*t (expected uniformised events) above which the
 * Tavaré series is used instead of uniformisation. */
extern double SX_UNIF_MAX_LAMBDA_T;

/* Poisson(lambda_t) weights over [m_lo, m_hi], normalised; stored in ws->w[m - m_lo]. */
void sx_poisson_weights(sx_unif_ws *ws, double lambda_t, int *m_lo, int *m_hi);

/* out[k] = sum_j v[j] P_{j,k}(t), k = 0..fmax (v over 0..fmax). N: haploid population size. */
void sx_unif_forward(sx_unif_ws *ws, int fmax, double t, double N, const double *v, double *out);

/* out[k - o] = P_{k,o}(t) for k in [o..kmax]. */
void sx_unif_column(sx_unif_ws *ws, int kmax, int o, double t, double N, double *out);

/* P_{i,o}(t) (single value). */
double sx_prob_coal(sx_unif_ws *ws, int i, int o, double t, double N);

/* Sample the i-o successive waiting times of a branch of length T conditioned
 * on going from i to o lineages.  waits[j] is the time from the previous event
 * (or the branch start) to the j-th coalescence.  Uses gsl_rng draws.
 * Returns SX_NO_ERROR or an error code if the fallback failed. */
long sx_sample_bounded_waits(sx_unif_ws *ws, int i, int o, double T, double N, float epsilon, gsl_rng *r, double *waits);

/* Tavaré alternating series in double with a cancellation check; falls back to
 * adaptive-precision MPFR.  Exposed for tests. */
double sx_tavare_prob(sx_unif_ws *ws, int i, int o, double t, double N);
/* Same, but returns 0 without escalating when |P| is provably below `neglig`. */
double sx_tavare_prob_bounded(sx_unif_ws *ws, int i, int o, double t, double N, double neglig);

#endif
