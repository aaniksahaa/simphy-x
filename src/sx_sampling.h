/* Sampling units ("f:1", "u:10,20", "ln:m,s", ...) — port of SimPhy's sampling.c
 * with identical semantics, including the dependency ordering used when several
 * units are sampled together (a stable sort by dependency depth). */
#ifndef SX_SAMPLING_H
#define SX_SAMPLING_H
#include "sx_common.h"
#include "sx_rng.h"

enum SX_DISTRIBUTIONS { SX_FIXED = 0, SX_UNIFORM = 1, SX_NORMAL = 2, SX_EXPONENTIAL = 3, SX_GAMMA = 4,
                        SX_LOGNORMAL = 5, SX_LOGUNIFORM = 6, SX_LOGNORMAL_MULT = 7 };
enum SX_TYPES { SX_UI = 0, SX_D = 1, SX_SU = 2 };

typedef struct sx_sunit sx_sunit;
typedef union { int i; double d; } sx_ui_d;
typedef union { int i; double d; sx_sunit *p; } sx_u_idsu;

struct sx_sunit {
    int distribution_code;
    sx_u_idsu params[5];
    int params_type[5];
    sx_ui_d value;
    int vtype;
    int dependent_index;
};

typedef struct { char name[3]; sx_sunit *p; } sx_sduple;
typedef struct { sx_sduple *table; int n_duples; } sx_stable;

#define sx_get(a) ((a).vtype == SX_UI ? (double)(a).value.i : (a).value.d)

void sx_set_uint(sx_sunit *v, int value);
void sx_set_double(sx_sunit *v, double value);
void sx_set_pointer_uint(sx_sunit *v, sx_sunit *dep);
void sx_set_pointer_double(sx_sunit *v, sx_sunit *dep);

/* Sample n units from r in SimPhy's order. */
long sx_sample_distr(gsl_rng *r, int n, sx_sunit **units);
long sx_parse_sampling(const char *p, sx_sunit *sample, const sx_stable tab);
void sx_print_sampling(const sx_sunit *sample, char *buffer, const sx_stable tab);
int sx_is_variable(sx_sunit value);
int sx_is_sampling_set(sx_sunit value);

#endif
