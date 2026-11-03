/* SPDX-License-Identifier: GPL-3.0-only */
#include "sx_sampling.h"
#include <ctype.h>
#include <strings.h>

static double casted_double(sx_u_idsu var, int kind)
{
    if (kind == SX_SU) return (double)sx_get(*var.p);
    else if (kind == SX_D) return var.d;
    else return (double)var.i;
}
static int casted_int(sx_u_idsu var, int kind)
{
    if (kind == SX_SU) return (int)sx_get(*var.p);
    else if (kind == SX_D) return (int)var.d;
    else return var.i;
}

#define set_prop(a, b) do { if ((a)->vtype == SX_UI) (a)->value.i = (int)(b); else (a)->value.d = (double)(b); } while (0)

void sx_set_uint(sx_sunit *v, int value)
{
    v->distribution_code = SX_FIXED; v->dependent_index = 0; v->value.i = value;
    v->params[0].i = value; v->params_type[0] = SX_UI; v->vtype = SX_UI;
}
void sx_set_double(sx_sunit *v, double value)
{
    v->distribution_code = SX_FIXED; v->dependent_index = 0; v->value.d = value;
    v->params[0].d = value; v->params_type[0] = SX_D; v->vtype = SX_D;
}
void sx_set_pointer_uint(sx_sunit *v, sx_sunit *dep)
{
    v->distribution_code = SX_FIXED; v->dependent_index = -1; v->value.i = 0;
    v->params[0].p = dep; v->params_type[0] = SX_SU; v->vtype = SX_UI;
}
void sx_set_pointer_double(sx_sunit *v, sx_sunit *dep)
{
    v->distribution_code = SX_FIXED; v->dependent_index = -1; v->value.d = 0;
    v->params[0].p = dep; v->params_type[0] = SX_SU; v->vtype = SX_D;
}

static int measure_recurrency(sx_sunit *su, int *it)
{
    int mdepth = 0, depth;
    ++*it;
    if (*it >= MAX_IT)
        sx_check(SX_SETTINGS_ERROR, "Recurrent sampling distribution parameters. Please, check your simulation parameters!\n");
    if (su->dependent_index == -1) {
        for (int i = 0; i < 5; ++i) {
            if (su->params_type[i] == SX_SU) {
                depth = measure_recurrency(su->params[i].p, it);
                if (depth > mdepth) mdepth = depth;
            }
        }
        su->dependent_index = mdepth + 1;
    }
    return su->dependent_index;
}

static int depth_of(sx_sunit *u)
{
    int rec = 0;
    if (u->dependent_index == -1) u->dependent_index = measure_recurrency(u, &rec);
    return u->dependent_index;
}

long sx_sample_distr(gsl_rng *r, int n, sx_sunit **in_units)
{
    sx_sunit **units = sx_malloc(sizeof(*units) * (size_t)n);
    memcpy(units, in_units, sizeof(*units) * (size_t)n);
    /* stable insertion sort by dependency depth (equivalent to glibc's stable qsort) */
    for (int i = 1; i < n; ++i) {
        sx_sunit *key = units[i];
        int kd = depth_of(key), j = i - 1;
        while (j >= 0 && depth_of(units[j]) > kd) { units[j + 1] = units[j]; --j; }
        units[j + 1] = key;
    }
    for (int i = 0; i < n; ++i) {
        sx_sunit *v = units[i];
        double a, b;
        switch (v->distribution_code) {
        case SX_FIXED:
            if (v->vtype == SX_D) v->value.d = casted_double(v->params[0], v->params_type[0]);
            else v->value.i = casted_int(v->params[0], v->params_type[0]);
            break;
        case SX_UNIFORM:
            set_prop(v, gsl_ran_flat(r, casted_double(v->params[0], v->params_type[0]),
                                     v->vtype == SX_UI ? casted_double(v->params[1], v->params_type[1]) + 1
                                                       : casted_double(v->params[1], v->params_type[1])));
            break;
        case SX_NORMAL:
            set_prop(v, (gsl_ran_gaussian_ziggurat(r, casted_double(v->params[1], v->params_type[1])) + casted_double(v->params[0], v->params_type[0])));
            break;
        case SX_EXPONENTIAL:
            set_prop(v, gsl_ran_exponential(r, 1 / casted_double(v->params[0], v->params_type[0])));
            break;
        case SX_GAMMA:
            set_prop(v, gsl_ran_gamma(r, casted_double(v->params[0], v->params_type[0]), casted_double(v->params[1], v->params_type[1])));
            break;
        case SX_LOGNORMAL:
            set_prop(v, gsl_ran_lognormal(r, casted_double(v->params[0], v->params_type[0]), casted_double(v->params[1], v->params_type[1])));
            break;
        case SX_LOGUNIFORM:
            a = log(casted_double(v->params[0], v->params_type[0]));
            b = log(casted_double(v->params[1], v->params_type[1]));
            set_prop(v, exp(a + (b - a) * gsl_rng_uniform(r)));
            break;
        case SX_LOGNORMAL_MULT:
            set_prop(v, (gsl_ran_lognormal(r, casted_double(v->params[0], v->params_type[0]), casted_double(v->params[1], v->params_type[1])) * casted_double(v->params[2], v->params_type[2])));
            break;
        default:
            free(units);
            return SX_UNEXPECTED_VALUE;
        }
    }
    free(units);
    return SX_NO_ERROR;
}

long sx_parse_sampling(const char *p, sx_sunit *sample, const sx_stable tab)
{
    char *buffer = sx_calloc(strlen(p) + 1, 1), *last_parsed = NULL;
    int in_parsing = 0, is_pointer = 0, b_index = 0, p_index = 0, i = 0, j = 0, p_found = 0, n_p = 0;
    long temp_int = 0;
    long rc = SX_NO_ERROR;

    switch (toupper((unsigned char)*p)) {
    case 'F': sample->distribution_code = SX_FIXED; n_p = 1; break;
    case 'U': sample->distribution_code = SX_UNIFORM; n_p = 2; break;
    case 'N': sample->distribution_code = SX_NORMAL; n_p = 2; break;
    case 'E': sample->distribution_code = SX_EXPONENTIAL; n_p = 1; break;
    case 'G': sample->distribution_code = SX_GAMMA; n_p = 2; break;
    case 'L':
        switch (toupper((unsigned char)*(p + 1))) {
        case 'N': sample->distribution_code = SX_LOGNORMAL; n_p = 2; ++p; break;
        case 'U': sample->distribution_code = SX_LOGUNIFORM; n_p = 2; ++p; break;
        default: free(buffer); return SX_SETTINGS_ERROR;
        }
        break;
    case 'S':
        switch (toupper((unsigned char)*(p + 1))) {
        case 'L': sample->distribution_code = SX_LOGNORMAL_MULT; n_p = 3; ++p; break;
        default: free(buffer); return SX_SETTINGS_ERROR;
        }
        break;
    default: free(buffer); return SX_SETTINGS_ERROR;
    }
    ++p;
    if (*p != ':') { free(buffer); return SX_SETTINGS_ERROR; }

    do {
        ++i; ++p;
        if (p_index > 4) { free(buffer); return SX_SETTINGS_ERROR; }
        switch (in_parsing) {
        case 0:
            in_parsing = 1; b_index = 1; buffer[0] = *p;
            if (!(isdigit((unsigned char)*p) || *p == '.' || *p == '-' || *p == '+')) is_pointer = 1;
            break;
        default:
            switch (*p) {
            case ',': case '\0': case '/': case '\n':
                in_parsing = 0; buffer[b_index] = '\0';
                switch (is_pointer) {
                default:
                    temp_int = strtol(buffer, &last_parsed, 10);
                    if ((temp_int >= INT_MIN && temp_int <= INT_MAX) && (*last_parsed == '\0' || *last_parsed == ',' || *last_parsed == '/' || *last_parsed == '\n')) {
                        sample->params[p_index].i = (int)temp_int;
                        sample->params_type[p_index] = SX_UI;
                    } else {
                        sample->params[p_index].d = strtod(buffer, &last_parsed);
                        if (!(*last_parsed == '\0' || *last_parsed == ',' || *last_parsed == '/' || *last_parsed == '\n')) { free(buffer); return SX_SETTINGS_ERROR; }
                        sample->params_type[p_index] = SX_D;
                    }
                    ++p_index;
                    break;
                case 1:
                    if (b_index > 2) { free(buffer); return SX_SETTINGS_ERROR; }
                    p_found = 0;
                    for (j = 0; j < tab.n_duples; ++j) {
                        if (strcasecmp(buffer, tab.table[j].name) == 0) {
                            sample->params[p_index].p = tab.table[j].p;
                            sample->params_type[p_index] = SX_SU;
                            sample->dependent_index = -1;
                            p_found = 1;
                            break;
                        }
                    }
                    if (p_found == 0) { free(buffer); return SX_SETTINGS_ERROR; }
                    ++p_index;
                    break;
                }
                is_pointer = 0;
                break;
            case ' ':
                break;
            default:
                buffer[b_index] = *p; ++b_index;
                break;
            }
            break;
        }
    } while (*p != '\0' && *p != '/' && i <= MAX_IT);
    if (i >= MAX_IT) rc = MAX_IT;
    else if (p_index != n_p || (*p == '/' && *(p + 1) != '/')) rc = SX_SETTINGS_ERROR;
    if (rc == SX_NO_ERROR && sample->distribution_code == SX_FIXED) {
        if (sample->params_type[0] == SX_UI) set_prop(sample, sample->params[0].i);
        else if (sample->params_type[0] == SX_D) set_prop(sample, sample->params[0].d);
    }
    free(buffer);
    return rc;
}

void sx_print_sampling(const sx_sunit *sample, char *buffer, const sx_stable tab)
{
    int n_param = 0, i, j;
    switch (sample->distribution_code) {
    case SX_FIXED: sprintf(buffer, "Fixed "); n_param = 1; break;
    case SX_UNIFORM: sprintf(buffer, "Uniform["); n_param = 2; break;
    case SX_NORMAL: sprintf(buffer, "Normal("); n_param = 2; break;
    case SX_EXPONENTIAL: sprintf(buffer, "Exp("); n_param = 1; break;
    case SX_GAMMA: sprintf(buffer, "Gamma("); n_param = 2; break;
    case SX_LOGNORMAL: sprintf(buffer, "LogN("); n_param = 2; break;
    case SX_LOGUNIFORM: sprintf(buffer, "LogU["); n_param = 2; break;
    case SX_LOGNORMAL_MULT: sprintf(buffer, "LogN("); n_param = 3; break;
    default: break;
    }
    for (i = 0; i < n_param; ++i) {
        switch (sample->params_type[i]) {
        case SX_UI: sprintf(buffer + strlen(buffer), "%u,", (unsigned)sample->params[i].i); break;
        case SX_D: sprintf(buffer + strlen(buffer), "%e,", sample->params[i].d); break;
        case SX_SU:
            for (j = 0; j <= tab.n_duples; ++j)
                if (tab.table[j].p == sample->params[i].p) break;
            sprintf(buffer + strlen(buffer), "%s,", tab.table[j].name);
            break;
        default: sx_check(SX_UNEXPECTED_VALUE, ""); break;
        }
    }
    switch (sample->distribution_code) {
    default: sprintf(buffer + strlen(buffer) - 1, ")"); break;
    case SX_UNIFORM:
        if (sample->vtype == SX_UI) sprintf(buffer + strlen(buffer) - 1, "]");
        else sprintf(buffer + strlen(buffer) - 1, ")");
        /* fallthrough (as in SimPhy) */
    case SX_FIXED: break;
    }
}

int sx_is_variable(sx_sunit value)
{
    switch (value.distribution_code) {
    case SX_FIXED:
        switch (value.params_type[0]) {
        case SX_SU: return sx_is_variable(*value.params[0].p);
        case SX_UI: return value.params[0].i == 0 ? 0 : 1;
        case SX_D: return value.params[0].d == 0 ? 0 : 1;
        }
        break;
    default: return 1;
    }
    return 0;
}

int sx_is_sampling_set(sx_sunit value)
{
    if (value.vtype == SX_UI)
        return !((value.distribution_code == SX_FIXED) && (value.value.i == 0) && (value.dependent_index == 0));
    return !((value.distribution_code == SX_FIXED) && (value.value.d == 0) && (value.dependent_index == 0));
}
