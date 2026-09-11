#include "sx_rng.h"
#include "sx_common.h"

gsl_rng *sx_rng_compat_new(unsigned long seed)
{
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    if (!r) sx_fatal("Cannot allocate the random number generator");
    gsl_rng_set(r, seed);
    return r;
}

/* ---- xoshiro256** ---- */
typedef struct { uint64_t s[4]; } xo_state;

static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static uint64_t splitmix64(uint64_t *x)
{
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline uint64_t xo_next(xo_state *st)
{
    uint64_t *s = st->s;
    const uint64_t result = rotl(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
}

static void xo_seed64(xo_state *st, uint64_t seed)
{
    uint64_t x = seed;
    for (int i = 0; i < 4; ++i) st->s[i] = splitmix64(&x);
    if (!(st->s[0] | st->s[1] | st->s[2] | st->s[3])) st->s[0] = 1;
}

/* gsl_rng_type interface: get() returns 32 bits (upper half of the 64-bit
 * output), get_double() uses 53 bits.  The seed passed by gsl_rng_set is only
 * unsigned long; we store the full 64 bits via sx_rng_fast_reseed. */
static void xo_set(void *vstate, unsigned long seed) { xo_seed64((xo_state *)vstate, (uint64_t)seed); }
static unsigned long xo_get(void *vstate) { return (unsigned long)(xo_next((xo_state *)vstate) >> 32); }
static double xo_get_double(void *vstate) { return (double)(xo_next((xo_state *)vstate) >> 11) * 0x1.0p-53; }

static const gsl_rng_type xo_type = {
    "xoshiro256**",
    0xffffffffUL,
    0,
    sizeof(xo_state),
    &xo_set,
    &xo_get,
    &xo_get_double
};
const gsl_rng_type *sx_rng_xoshiro = &xo_type;

uint64_t sx_stream_seed(uint64_t seed, uint64_t replicate, uint64_t locus, uint64_t tag)
{
    uint64_t x = seed ^ 0xD1B54A32D192ED03ULL;
    uint64_t h = splitmix64(&x);
    x ^= replicate * 0x9E3779B97F4A7C15ULL; h ^= splitmix64(&x);
    x ^= locus * 0xC2B2AE3D27D4EB4FULL; h ^= splitmix64(&x);
    x ^= tag * 0x165667B19E3779F9ULL; h ^= splitmix64(&x);
    return h;
}

gsl_rng *sx_rng_fast_new(uint64_t stream_seed)
{
    gsl_rng *r = gsl_rng_alloc(sx_rng_xoshiro);
    if (!r) sx_fatal("Cannot allocate the random number generator");
    sx_rng_fast_reseed(r, stream_seed);
    return r;
}

void sx_rng_fast_reseed(gsl_rng *r, uint64_t stream_seed)
{
    xo_seed64((xo_state *)gsl_rng_state(r), stream_seed);
}
