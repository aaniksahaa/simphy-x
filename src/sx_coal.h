/* Compat numerics: exact port of SimPhy's bounded multilocus coalescent
 * machinery (Kahan log-scale Tavaré series, modified Brent root finding, and
 * the 512-bit MPFR fallback).  Every floating-point expression is kept in
 * SimPhy's form so that results are bit-identical. */
#ifndef SX_COAL_H
#define SX_COAL_H
#include "sx_common.h"
#include "sx_rng.h"

double sx_logscale_fact(int n);
double sx_biomcmc_log1p(double x);
double sx_logscale_add(double logx, double logy);
double sx_kahan_logscale_add(double logx, double logy, double *comp_sum);
double sx_logscale_sub(double logx, double logy);
void sx_skahan_sum(double *sum, double input, double *compensation);
double sx_vkahan_sum(const double *array, int n);
size_t sx_sample_ndoubles(size_t n, const double *array, gsl_rng *seed);

double sx_kahan_logscale_prob_coal(int i_lin, int o_lin, double bound_time, int Ne);

/* Double-precision attempt (SampleCoalTimeMLCFromXtoYLineages). */
double sx_sample_coal_time_mlc(int i_lin, int o_lin, double bound_time, int pop_size, double brent_epsilon, double density, int verbosity);

/* MPFR chain: SimPhy redraws all n_coals waiting times with MPFR when the
 * double attempt failed.  u[] are the uniforms in draw order; out[] receives
 * the sampled times.  Deterministic, thread-safe (own MPFR variables). */
typedef struct {
    int n_coals, n_ilin, n_olin, pop_size;
    double n_gen, min_ngen;
    float epsilon;
    int verbosity;
    double *u;      /* n_coals uniforms */
    double *out;    /* n_coals results */
} sx_mpfr_chain_job;
void sx_mpfr_chain_run(sx_mpfr_chain_job *job);
/* 0: abort like SimPhy 1.0.2 when the 512-bit CDF leaves [0,1]; 1: retry with doubled precision. */
extern int sx_mpfr_rescue;
extern unsigned long sx_mpfr_rescue_count;
double sx_mpfr_sample_coal_time_mlc(int i_lin, int o_lin, double bound_time, int pop_size, double brent_epsilon, double density, int verbosity, int precision);
/* Diagnostics: SimPhy's conditional CDF of the next coalescence time evaluated at w (MPFR / double versions). */
double sx_debug_mpfr_cdf(int i_lin, int o_lin, double bound_time, int pop_size, double w, int precision);
double sx_debug_double_cdf(int i_lin, int o_lin, double bound_time, int pop_size, double w);

#endif
