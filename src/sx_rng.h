/* SPDX-License-Identifier: GPL-3.0-only */
#ifndef SX_RNG_H
#define SX_RNG_H
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdint.h>

/* Compat stream: SimPhy's generator (gsl_rng_default after gsl_rng_env_setup,
 * i.e. mt19937 unless GSL_RNG_TYPE says otherwise), seeded with gsl_rng_set. */
gsl_rng *sx_rng_compat_new(unsigned long seed);

/* Fast-mode streams: xoshiro256** exposed as a gsl_rng_type, seeded from a
 * 64-bit value through splitmix64.  Stream ids are derived deterministically
 * from (seed, replicate, locus). */
extern const gsl_rng_type *sx_rng_xoshiro;
uint64_t sx_stream_seed(uint64_t seed, uint64_t replicate, uint64_t locus, uint64_t tag);
gsl_rng *sx_rng_fast_new(uint64_t stream_seed);
void sx_rng_fast_reseed(gsl_rng *r, uint64_t stream_seed);

#endif
