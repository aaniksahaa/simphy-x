#include "sx_unif.h"
#include "sx_coal.h"
#include <mpfr.h>

double SX_UNIF_MAX_LAMBDA_T = 20000.0;

void sx_unif_ws_init(sx_unif_ws *ws) { memset(ws, 0, sizeof(*ws)); }
void sx_unif_ws_free(sx_unif_ws *ws) { free(ws->w); free(ws->b); free(ws->v0); free(ws->v1); free(ws->tau); memset(ws, 0, sizeof(*ws)); }

static inline double lam(int k, double N) { return (double)k * (double)(k - 1) / (2.0 * N); }

static void ws_w(sx_unif_ws *ws, int n) { if (n > ws->w_cap) { ws->w = sx_realloc(ws->w, sizeof(double) * (size_t)n); ws->w_cap = n; } }
static void ws_v(sx_unif_ws *ws, int n) { if (n > ws->v_cap) { ws->v0 = sx_realloc(ws->v0, sizeof(double) * (size_t)n); ws->v1 = sx_realloc(ws->v1, sizeof(double) * (size_t)n); ws->v_cap = n; } }
static void ws_b(sx_unif_ws *ws, size_t n) { if (n > ws->b_cap) { ws->b = sx_realloc(ws->b, sizeof(double) * n); ws->b_cap = n; } }
static void ws_tau(sx_unif_ws *ws, int n) { if (n > ws->tau_cap) { ws->tau = sx_realloc(ws->tau, sizeof(double) * (size_t)n); ws->tau_cap = n; } }

/* Upper truncation point of the event count: the chain needs at least J
 * jumps, and the Poisson tail beyond mode + 12 sqrt + 50 is negligible.  The
 * callers verify that the last contribution is negligible and extend if not. */
static int poisson_hi(double lt, int J) { return (int)(lt + 12.0 * sqrt(lt + (double)J) + 50.0) + J; }

/* log Poisson(m; lt) for all m in [0, hi], stored in ws->w. */
static void log_poisson(sx_unif_ws *ws, double lt, int hi)
{
    ws_w(ws, hi + 1);
    if (lt <= 0) { for (int m = 0; m <= hi; ++m) ws->w[m] = (m == 0) ? 0.0 : -INFINITY; return; }
    double llt = log(lt), lf = 0;
    for (int m = 0; m <= hi; ++m) {
        if (m > 0) lf += log((double)m);
        ws->w[m] = m * llt - lt - lf;
    }
}

void sx_poisson_weights(sx_unif_ws *ws, double lt, int *m_lo, int *m_hi)
{
    int hi = lt <= 0 ? 0 : poisson_hi(lt, 0);
    log_poisson(ws, lt, hi);
    *m_lo = 0; *m_hi = hi;
}

#define SX_TAIL_LOG 46.0   /* stop once contributions fell 1e-20 below the largest one */

/* Online log-scaled accumulation: acc holds values relative to exponent E. */
typedef struct { double E; int init; } lacc;

static inline void lacc_add(lacc *a, double *acc, int n, double e_m, const double *vec)
{
    if (!a->init) { a->E = e_m; a->init = 1; for (int k = 0; k < n; ++k) acc[k] = vec[k]; return; }
    if (e_m > a->E + 200.0) {
        double f = exp(a->E - e_m);
        for (int k = 0; k < n; ++k) acc[k] = acc[k] * f + vec[k];
        a->E = e_m;
        return;
    }
    double f = exp(e_m - a->E);
    if (f == 0.0) return;
    for (int k = 0; k < n; ++k) acc[k] += f * vec[k];
}

/* Rescale a running vector so its max is 1 when it drifted below 1e-100 (or above 1); returns log of the factor removed. */
static inline double renorm(double *v, int n)
{
    double mx = 0;
    for (int k = 0; k < n; ++k) if (v[k] > mx) mx = v[k];
    if (mx <= 0 || (mx >= 1e-100 && mx <= 1.0)) return 0.0;
    double l = log(mx), f = 1.0 / mx;
    for (int k = 0; k < n; ++k) v[k] *= f;
    return l;
}

/* ---- forward transport: out = v e^{Qt} over states 0..fmax ---- */
static void unif_forward_core(sx_unif_ws *ws, int fmax, double t, double N, const double *v, double *out)
{
    double Lam = lam(fmax, N), lt = Lam * t;
    int n = fmax + 1, J = fmax;
    int hi = poisson_hi(lt, J);
    for (;;) {
        log_poisson(ws, lt, hi);
        ws_v(ws, n);
        double *cur = ws->v0, *nxt = ws->v1;
        memcpy(cur, v, sizeof(double) * (size_t)n);
        double L = renorm(cur, n);
        double *best = sx_malloc(sizeof(double) * (size_t)n); for (int k = 0; k < n; ++k) best[k] = -INFINITY;
        int converged = 0;
        lacc a = { 0, 0 };
        for (int k = 0; k < n; ++k) out[k] = 0;
        for (int m = 0; m <= hi; ++m) {
            double e = ws->w[m] + L;
            lacc_add(&a, out, n, e, cur);
            int all_small = 1;
            for (int k = 0; k < n; ++k) {
                double lc = cur[k] > 0 ? e + log(cur[k]) : -INFINITY;
                if (lc > best[k]) best[k] = lc;
                if (best[k] > -INFINITY && lc >= best[k] - SX_TAIL_LOG) all_small = 0;
            }
            if (m >= J && m > lt && all_small) { converged = 1; hi = m; break; }
            if (m == hi) break;
            for (int k = 0; k <= fmax; ++k) {
                double stay = 1.0 - lam(k, N) / Lam;
                double in = (k + 1 <= fmax) ? cur[k + 1] * (lam(k + 1, N) / Lam) : 0.0;
                nxt[k] = cur[k] * stay + in;
            }
            double *tmp = cur; cur = nxt; nxt = tmp;
            L += renorm(cur, n);
        }
        free(best);
        if (converged) {
            double f = a.init ? exp(a.E) : 0.0;
            for (int k = 0; k < n; ++k) out[k] *= f;
            return;
        }
        hi = hi * 2 + 100;   /* tail not yet negligible: extend and redo */
    }
}

/* Adaptive uniformisation: the pure-death chain leaves high states quickly,
 * so the active state space and the uniformisation rate are shrunk piecewise.
 * Total cost is O(fmax^2) regardless of the branch length. */
#define SX_ADAPT_BUDGET 4096.0    /* max uniformised events per piece */
#define SX_ADAPT_TRUNC 1e-60      /* relative mass below which high states are dropped */

void sx_unif_forward(sx_unif_ws *ws, int fmax, double t, double N, const double *v, double *out)
{
    if (fmax <= 1 || t <= 0) { memcpy(out, v, sizeof(double) * (size_t)(fmax + 1)); return; }
    double *cur = sx_malloc(sizeof(double) * (size_t)(fmax + 1));
    double *nxt = sx_malloc(sizeof(double) * (size_t)(fmax + 1));
    memcpy(cur, v, sizeof(double) * (size_t)(fmax + 1));
    double t_rem = t;
    int kmax = fmax;
    int pieces = 0;
    while (t_rem > 0) {
        /* active range */
        double mx = 0;
        for (int k = 0; k <= kmax; ++k) if (cur[k] > mx) mx = cur[k];
        if (mx <= 0) break;
        while (kmax > 1 && cur[kmax] <= SX_ADAPT_TRUNC * mx) { cur[kmax] = 0; --kmax; }
        if (kmax <= 1) break;                       /* nothing left to coalesce */
        double Lam = lam(kmax, N);
        double budget = 2.0 * kmax; if (budget < 128.0) budget = 128.0; if (budget > SX_ADAPT_BUDGET) budget = SX_ADAPT_BUDGET;
        double delta = budget / Lam;
        if (delta >= t_rem) delta = t_rem;
        unif_forward_core(ws, kmax, delta, N, cur, nxt);
        for (int k = 0; k <= kmax; ++k) cur[k] = nxt[k];
        t_rem -= delta;
        ++pieces;
    }
    memcpy(out, cur, sizeof(double) * (size_t)(fmax + 1));
    ws->n_forward_tavare += (unsigned long)pieces;   /* reused as "pieces" counter */
    free(cur); free(nxt);
}

/* ---- column: out[k-o] = P_{k,o}(t), k in [o..kmax] ---- */
static void unif_column_core(sx_unif_ws *ws, int kmax, int o, double t, double N, double *out)
{
    int n = kmax - o + 1, J = n - 1;
    double Lam = lam(kmax, N), lt = Lam * t;
    int hi = poisson_hi(lt, J);
    for (;;) {
        log_poisson(ws, lt, hi);
        ws_v(ws, n);
        double *cur = ws->v0, *nxt = ws->v1;
        for (int k = 0; k < n; ++k) { cur[k] = 0; out[k] = 0; }
        cur[0] = 1.0;
        double L = 0;
        double *best = sx_malloc(sizeof(double) * (size_t)n); for (int k = 0; k < n; ++k) best[k] = -INFINITY;
        int converged = 0;
        lacc a = { 0, 0 };
        for (int m = 0; m <= hi; ++m) {
            double e = ws->w[m] + L;
            lacc_add(&a, out, n, e, cur);
            int all_small = 1;
            for (int k = 0; k < n; ++k) {
                double lc = cur[k] > 0 ? e + log(cur[k]) : -INFINITY;
                if (lc > best[k]) best[k] = lc;
                if (best[k] > -INFINITY && lc >= best[k] - SX_TAIL_LOG) all_small = 0;
            }
            if (m >= J && m > lt && all_small) { converged = 1; hi = m; break; }
            if (m == hi) break;
            for (int k = 0; k < n; ++k) {
                double lk = lam(o + k, N) / Lam;
                nxt[k] = cur[k] * (1.0 - lk) + (k > 0 ? cur[k - 1] * lk : 0.0);
            }
            double *tmp = cur; cur = nxt; nxt = tmp;
            L += renorm(cur, n);
        }
        free(best);
        if (converged) {
            double f = a.init ? exp(a.E) : 0.0;
            for (int k = 0; k < n; ++k) out[k] *= f;
            return;
        }
        hi = hi * 2 + 100;
    }
}

void sx_unif_column(sx_unif_ws *ws, int kmax, int o, double t, double N, double *out)
{
    int n = kmax - o + 1;
    if (n <= 0) return;
    double Lam = lam(kmax, N);
    if (t <= 0) { for (int k = 0; k < n; ++k) out[k] = (k == 0) ? 1.0 : 0.0; return; }
    if (Lam * t <= SX_UNIF_MAX_LAMBDA_T) { unif_column_core(ws, kmax, o, t, N, out); return; }
    for (int k = 0; k < n; ++k) out[k] = sx_tavare_prob(ws, o + k, o, t, N);
}

double sx_prob_coal(sx_unif_ws *ws, int i, int o, double t, double N)
{
    if (i < o || o < 1) return 0.0;
    if (i == o) return exp(-lam(o, N) * t);
    if (t <= 0) return 0.0;
    double Lam = lam(i, N);
    if (Lam * t <= SX_UNIF_MAX_LAMBDA_T) {
        double *tmp = sx_malloc(sizeof(double) * (size_t)(i - o + 1));
        unif_column_core(ws, i, o, t, N, tmp);
        double p = tmp[i - o];
        free(tmp);
        return p;
    }
    return sx_tavare_prob(ws, i, o, t, N);
}

/* ---- Tavaré series ----
 * P_{i,o}(t) = sum_{k=o}^{i} exp(-lam_k t) c_k (-1)^{k-o},
 * c_o = prod_{j<o} (o+j)(i-j)/(i+j) / o!,
 * c_k / c_{k-1} = (2k-1)(o+k-2)(i-k+1) / ((2k-3)(k-o)(i+k-1)). */
static double tavare_double(int i, int o, double t, double N, double *ratio)
{
    double pos = 0, neg = 0, cpos = 0, cneg = 0;
    double logc = 0;
    for (int j = 0; j < o; ++j) logc += log((double)(o + j)) + log((double)(i - j)) - log((double)(i + j));
    logc -= sx_logscale_fact(o);
    double term = exp(-lam(o, N) * t + logc);
    pos = term;
    double maxmag = term;
    for (int k = o + 1; k <= i; ++k) {
        logc += log((2.0 * k - 1) / (2.0 * k - 3)) + log((double)(o + k - 2)) + log((double)(i - k + 1)) - log((double)(k - o)) - log((double)(i + k - 1));
        term = exp(-lam(k, N) * t + logc);
        if (term > maxmag) maxmag = term;
        if ((k - o) % 2 == 0) sx_skahan_sum(&pos, term, &cpos); else sx_skahan_sum(&neg, term, &cneg);
    }
    double res = pos - neg;
    *ratio = (res > 0) ? maxmag / res : INFINITY;
    return res;
}

static double tavare_mpfr(int i, int o, double t, double N, int prec)
{
    mpfr_t sum, term, c, e;
    mpfr_inits2(prec, sum, term, c, e, (mpfr_ptr)0);
    mpfr_set_zero(sum, 1);
    mpfr_set_si(c, 1, MPFR_RNDN);
    for (int j = 0; j < o; ++j) {
        mpfr_mul_si(c, c, (long)(o + j) * (long)(i - j), MPFR_RNDN);
        mpfr_div_si(c, c, i + j, MPFR_RNDN);
    }
    mpfr_fac_ui(e, (unsigned long)o, MPFR_RNDN);
    mpfr_div(c, c, e, MPFR_RNDN);
    for (int k = o; k <= i; ++k) {
        if (k > o) {
            mpfr_mul_si(c, c, (long)(2 * k - 1) * (long)(o + k - 2), MPFR_RNDN);
            mpfr_mul_si(c, c, i - k + 1, MPFR_RNDN);
            mpfr_div_si(c, c, (long)(2 * k - 3) * (long)(k - o), MPFR_RNDN);
            mpfr_div_si(c, c, i + k - 1, MPFR_RNDN);
        }
        mpfr_set_si(e, -(long)k * (long)(k - 1), MPFR_RNDN);
        mpfr_mul_d(e, e, t, MPFR_RNDN);
        mpfr_div_d(e, e, 2.0 * N, MPFR_RNDN);
        mpfr_exp(e, e, MPFR_RNDN);
        mpfr_mul(term, e, c, MPFR_RNDN);
        if ((k - o) % 2 == 0) mpfr_add(sum, sum, term, MPFR_RNDN); else mpfr_sub(sum, sum, term, MPFR_RNDN);
    }
    double r = mpfr_get_d(sum, MPFR_RNDN);
    mpfr_clears(sum, term, c, e, (mpfr_ptr)0);
    return r;
}

double sx_tavare_prob_bounded(sx_unif_ws *ws, int i, int o, double t, double N, double neglig)
{
    if (i < o || o < 1) return 0.0;
    if (i == o) return exp(-lam(o, N) * t);
    double ratio;
    double p = tavare_double(i, o, t, N, &ratio);
    if (ratio <= 1e3 && p >= 0) return p;
    /* |P| <= sum of term magnitudes = ratio-free bound: maxmag * (i-o+1) is a crude bound */
    if (neglig > 0) {
        double bound = (p > 0 ? p : 0) ;
        /* the largest term magnitude bounds each of the (i-o+1) terms */
        double maxmag = (ratio < INFINITY && p > 0) ? ratio * p : INFINITY;
        if (maxmag < INFINITY) bound = maxmag * (i - o + 1);
        else {
            /* recompute the bound directly */
            double logc = 0, mm = 0;
            for (int j = 0; j < o; ++j) logc += log((double)(o + j)) + log((double)(i - j)) - log((double)(i + j));
            logc -= sx_logscale_fact(o);
            mm = exp(-lam(o, N) * t + logc);
            for (int k = o + 1; k <= i; ++k) {
                logc += log((2.0 * k - 1) / (2.0 * k - 3)) + log((double)(o + k - 2)) + log((double)(i - k + 1)) - log((double)(k - o)) - log((double)(i + k - 1));
                double term = exp(-lam(k, N) * t + logc);
                if (term > mm) mm = term;
            }
            bound = mm * (i - o + 1);
        }
        if (bound < neglig) return 0.0;
    }
    ws->n_mpfr++;
    double lr = (ratio > 1 && ratio < INFINITY) ? log2(ratio) : 200.0;
    int prec = 128;
    while (prec < lr + 64) prec *= 2;
    for (;;) {
        double a = tavare_mpfr(i, o, t, N, prec);
        double b = tavare_mpfr(i, o, t, N, prec * 2);
        if (b == 0 || fabs(a - b) <= 1e-13 * fabs(b) || prec >= 32768) return b;
        prec *= 4;
    }
}

double sx_tavare_prob(sx_unif_ws *ws, int i, int o, double t, double N)
{
    return sx_tavare_prob_bounded(ws, i, o, t, N, 0.0);
}

/* ---- bridge sampler ---- */
static int cmp_dbl(const void *a, const void *b) { double x = *(const double *)a, y = *(const double *)b; return (x < y) ? -1 : (x > y); }

static long fallback_inversion(sx_unif_ws *ws, int i, int o, double T, double N, float epsilon, gsl_rng *r, double *waits, int mpfr_only)
{
    ws->n_bridge_fallback++;
    int n_ev = i - o;
    double current = T; int avail = i;
    for (int j = 0; j < n_ev; ++j) {
        double w = mpfr_only ? -1.0 : sx_sample_coal_time_mlc(avail, o, current, (int)N, epsilon, gsl_rng_uniform_pos(r), 0);
        if (!(w > 0) || w >= current) {
            for (int jj = j; jj < n_ev; ++jj) {
                w = sx_mpfr_sample_coal_time_mlc(avail, o, current, (int)N, epsilon, gsl_rng_uniform_pos(r), 0, 512);
                waits[jj] = w; current -= w; --avail;
            }
            mpfr_free_cache();
            return SX_NO_ERROR;
        }
        waits[j] = w; current -= w; --avail;
    }
    return SX_NO_ERROR;
}

/* memory bound for the stored backward vectors of a direct bridge */
#define SX_BRIDGE_MAX_DOUBLES 1e6

static int bridge_direct_ok(int i, int o, double T, double N)
{
    int n = i - o + 1, J = n - 1;
    double LT = lam(i, N) * T;
    if (LT > SX_UNIF_MAX_LAMBDA_T) return 0;
    return (double)(poisson_hi(LT, J) + 1) * (double)n <= SX_BRIDGE_MAX_DOUBLES;
}

static long bridge_direct(sx_unif_ws *ws, int i, int o, double T, double N, float epsilon, gsl_rng *r, double *waits);

/* Divide and conquer for many expected events: sample X(t1) for a short first
 * piece, bridge the first piece directly, recurse on the remainder. */
static long bridge_dc(sx_unif_ws *ws, int i, int o, double T, double N, float epsilon, gsl_rng *r, double *waits)
{
    int n_ev = i - o, filled = 0;
    double carry = 0;         /* time since the last coalescence */
    int cur_i = i; double t_rem = T;
    int level = 0;
    while (cur_i > o) {
        if (bridge_direct_ok(cur_i, o, t_rem, N)) {
            long rc = bridge_direct(ws, cur_i, o, t_rem, N, epsilon, r, waits + filled);
            if (rc != SX_NO_ERROR) return rc;
            waits[filled] += carry;
            return SX_NO_ERROR;
        }
        if (++level > 200) return fallback_inversion(ws, cur_i, o, t_rem, N, epsilon, r, waits + filled, 0);
        int n = cur_i - o + 1, J = n - 1;
        double Lam = lam(cur_i, N);
        /* first-piece event budget from the memory bound and the series threshold */
        double lt1 = SX_BRIDGE_MAX_DOUBLES / n - 12.0 * sqrt(SX_BRIDGE_MAX_DOUBLES / n) - 60.0 - J;
        if (lt1 > SX_UNIF_MAX_LAMBDA_T) lt1 = SX_UNIF_MAX_LAMBDA_T;
        if (lt1 < 64.0) return fallback_inversion(ws, cur_i, o, t_rem, N, epsilon, r, waits + filled, 0);
        double t1 = lt1 / Lam;
        if (t1 >= t_rem) t1 = t_rem * 0.5;
        /* row: P_{cur_i,k}(t1) for all k */
        double *row = sx_calloc((size_t)cur_i + 1, sizeof(double)), *e = sx_calloc((size_t)cur_i + 1, sizeof(double));
        e[cur_i] = 1.0;
        sx_unif_forward(ws, cur_i, t1, N, e, row);
        double rmax = 0; for (int k = o; k <= cur_i; ++k) if (row[k] > rmax) rmax = row[k];
        /* band-limited column entries G(k) = P_{k,o}(t_rem - t1) via adaptive rows */
        double *w = sx_calloc((size_t)cur_i + 1, sizeof(double));
        double s2 = t_rem - t1, wmax = 0;
        for (int k = o; k <= cur_i; ++k) {
            if (row[k] <= 1e-30 * rmax) continue;
            double g;
            if (k == o) g = exp(-lam(o, N) * s2);
            else if (lam(k, N) * s2 <= SX_UNIF_MAX_LAMBDA_T) g = sx_prob_coal(ws, k, o, s2, N);
            else {
                double *ek = sx_calloc((size_t)k + 1, sizeof(double)), *rk = sx_calloc((size_t)k + 1, sizeof(double));
                ek[k] = 1.0; sx_unif_forward(ws, k, s2, N, ek, rk); g = rk[o]; free(ek); free(rk);
            }
            w[k] = row[k] * g;
            if (w[k] > wmax) wmax = w[k];
        }
        free(row); free(e);
        if (!(wmax > 0)) { free(w); return fallback_inversion(ws, cur_i, o, t_rem, N, epsilon, r, waits + filled, 0); }
        double tot = 0; for (int k = o; k <= cur_i; ++k) tot += w[k];
        double u = gsl_rng_uniform(r) * tot, acc = 0; int ksel = cur_i;
        for (int k = o; k <= cur_i; ++k) { acc += w[k]; if (acc >= u) { ksel = k; break; } }
        free(w);
        /* first piece: cur_i -> ksel over t1 */
        int ev1 = cur_i - ksel;
        if (ev1 > 0) {
            long rc = bridge_direct(ws, cur_i, ksel, t1, N, epsilon, r, waits + filled);
            if (rc != SX_NO_ERROR) return rc;
            waits[filled] += carry;
            double used = 0; for (int j = 0; j < ev1; ++j) used += waits[filled + j];
            used -= carry;
            carry = t1 - used;
            filled += ev1;
        } else carry += t1;
        cur_i = ksel; t_rem -= t1;
    }
    (void)n_ev;
    return SX_NO_ERROR;
}

long sx_sample_bounded_waits(sx_unif_ws *ws, int i, int o, double T, double N, float epsilon, gsl_rng *r, double *waits)
{
    if (i - o <= 0) return SX_NO_ERROR;
    if (bridge_direct_ok(i, o, T, N)) return bridge_direct(ws, i, o, T, N, epsilon, r, waits);
    return bridge_dc(ws, i, o, T, N, epsilon, r, waits);
}

static long bridge_direct(sx_unif_ws *ws, int i, int o, double T, double N, float epsilon, gsl_rng *r, double *waits)
{
    int n_ev = i - o;
    if (n_ev <= 0) return SX_NO_ERROR;
    double Lam = lam(i, N);
    double LT = Lam * T;
    int n = i - o + 1, J = n - 1;
    int hi = poisson_hi(LT, J);
    double *Lb = NULL, *b = NULL;
    int m_hi = 0;
    for (;;) {
        if ((double)(hi + 1) * (double)n > SX_BRIDGE_MAX_DOUBLES) { free(Lb); return fallback_inversion(ws, i, o, T, N, epsilon, r, waits, 0); }
        log_poisson(ws, LT, hi);
        ws_b(ws, (size_t)(hi + 1) * (size_t)n);
        b = ws->b;
        Lb = sx_realloc(Lb, sizeof(double) * (size_t)(hi + 1));
        for (int k = 0; k < n; ++k) b[k] = (k == 0) ? 1.0 : 0.0;
        Lb[0] = 0;
        double best = -INFINITY, lastc = -INFINITY;
        m_hi = hi;
        for (int m = 0; m <= hi; ++m) {
            if (m > 0) {
                const double *prev = b + (size_t)(m - 1) * n;
                double *cur = b + (size_t)m * n;
                for (int k = 0; k < n; ++k) {
                    double lk = lam(o + k, N) / Lam;
                    cur[k] = prev[k] * (1.0 - lk) + (k > 0 ? prev[k - 1] * lk : 0.0);
                }
                Lb[m] = Lb[m - 1] + renorm(cur, n);
            }
            double bm = b[(size_t)m * n + (n - 1)];
            lastc = bm > 0 ? ws->w[m] + Lb[m] + log(bm) : -INFINITY;
            ws->w[m] = lastc;                   /* log joint weight of m events */
            if (lastc > best) best = lastc;
            if (m >= J && m > LT && lastc < best - SX_TAIL_LOG) { m_hi = m; break; }
        }
        if (lastc < best - SX_TAIL_LOG || best == -INFINITY) break;
        hi = hi * 2 + 100;
    }
    if (!(ws->w[0] > -INFINITY || m_hi > 0)) { free(Lb); return fallback_inversion(ws, i, o, T, N, epsilon, r, waits, 1); }
    ws->n_bridge++;
    double lmax = -INFINITY;
    for (int m = 0; m <= m_hi; ++m) if (ws->w[m] > lmax) lmax = ws->w[m];
    if (!(lmax > -INFINITY)) { free(Lb); return fallback_inversion(ws, i, o, T, N, epsilon, r, waits, 1); }
    double tot = 0;
    for (int m = 0; m <= m_hi; ++m) { double p = exp(ws->w[m] - lmax); ws->w[m] = p; tot += p; }
    double u = gsl_rng_uniform(r) * tot, acc = 0;
    int m_sel = m_hi;
    for (int m = 0; m <= m_hi; ++m) { acc += ws->w[m]; if (acc >= u) { m_sel = m; break; } }

    ws_tau(ws, m_sel + 1);
    double *tau = ws->tau;
    for (int e = 0; e < m_sel; ++e) tau[e] = gsl_rng_uniform_pos(r) * T;
    qsort(tau, (size_t)m_sel, sizeof(double), cmp_dbl);
    int x = i, jumps = 0;
    double last = 0;
    for (int e = 1; e <= m_sel; ++e) {
        int rem = m_sel - e;
        double b_cur = b[(size_t)(rem + 1) * n + (x - o)];
        if (x > o && b_cur > 0) {
            double p_jump = (lam(x, N) / Lam) * exp(Lb[rem] - Lb[rem + 1]) * b[(size_t)rem * n + (x - 1 - o)] / b_cur;
            if (gsl_rng_uniform(r) < p_jump) {
                waits[jumps++] = tau[e - 1] - last;
                last = tau[e - 1];
                --x;
            }
        }
    }
    free(Lb);
    if (ws->b_cap > (size_t)2.5e5) { free(ws->b); ws->b = NULL; ws->b_cap = 0; }   /* do not pin large scratch per thread */
    if (jumps != n_ev) return SX_UNEXPECTED_VALUE;
    return SX_NO_ERROR;
}
