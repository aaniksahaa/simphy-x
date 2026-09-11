/* glibc-compatible stable top-down mergesort (the algorithm qsort() uses on
 * glibc, msort_with_tmp).  Used instead of qsort so that results with
 * SimPhy's epsilon comparators are identical on every platform. */
#ifndef SX_MSORT_H
#define SX_MSORT_H
#include <stddef.h>
typedef int (*sx_cmp_fn)(const void *, const void *);
void sx_msort(void *base, size_t n, size_t size, sx_cmp_fn cmp);
/* SimPhy's Compare_DBL: equal when |a-b| < DBL_EPSILON. */
int sx_compare_dbl(const void *a, const void *b);
#endif
