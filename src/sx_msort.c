/* SPDX-License-Identifier: GPL-3.0-only */
#include "sx_msort.h"
#include "sx_common.h"

static void msort_rec(char *b, size_t n, size_t s, sx_cmp_fn cmp, char *tmp)
{
    if (n <= 1) return;
    size_t n1 = n / 2, n2 = n - n1;
    char *b1 = b, *b2 = b + n1 * s;
    msort_rec(b1, n1, s, cmp, tmp);
    msort_rec(b2, n2, s, cmp, tmp);
    char *t = tmp;
    while (n1 > 0 && n2 > 0) {
        if (cmp(b1, b2) <= 0) { memcpy(t, b1, s); b1 += s; --n1; }
        else { memcpy(t, b2, s); b2 += s; --n2; }
        t += s;
    }
    if (n1 > 0) memcpy(t, b1, n1 * s);
    memcpy(b, tmp, (n - n2) * s);
}

void sx_msort(void *base, size_t n, size_t size, sx_cmp_fn cmp)
{
    if (n <= 1) return;
    char *tmp = sx_malloc(n * size);
    msort_rec((char *)base, n, size, cmp, tmp);
    free(tmp);
}

int sx_compare_dbl(const void *n1, const void *n2)
{
    double a = *(const double *)n1, b = *(const double *)n2, dif;
    if (a - b > 0) dif = a - b; else dif = b - a;
    if (dif < DBL_EPSILON) return 0;
    return (a > b) ? 1 : -1;
}
