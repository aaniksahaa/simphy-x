/* simphy-x — common definitions shared by every module. */
#ifndef SX_COMMON_H
#define SX_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <stdarg.h>

/* Error codes: identical numeric values to SimPhy's enum ERRORS so that the
 * abort/continue semantics of ErrorReporter are preserved. */
enum {
    SX_MEM_ERROR = 0,
    SX_NO_ERROR = 1,
    SX_SETTINGS_ERROR = 2,
    SX_LOOP_ERROR = 3,
    SX_UNEXPECTED_VALUE = 4,
    SX_DB_ERROR = 5,
    SX_TERMINATE_NOERROR = 6,
    SX_IO_ERROR = -1
};

/* Locus-tree node kinds (SimPhy NODE_VALUES). */
enum { SP = 0, DUP = 1, LOSS = 2, TRFR = 3, RTRFR = 4, GC = 5, RGC = 6 };
/* Length units (SimPhy UNITS). */
enum { CU = 0, BL = 1, GL = 2, TL = 3 };

/* Run modes. */
enum { SX_MODE_COMPAT = 0, SX_MODE_FAST = 1 };

/* Tunable limits (SimPhy globals; overridable through the same env vars). */
extern int MAX_IT;
extern int MAX_NAME;
extern int MAX_LEAVES;
extern int NUM_BUFFER;
extern int IO_BUFFER;
extern char TEST_CHAR;
extern double FLOAT_PRECISION;

void sx_limits_from_env(void);

/* ErrorReporter equivalent: aborts on error codes, ignores NO_ERROR and
 * "pointer" values (anything outside the enum). */
void sx_check(long code, const char *msg);
void sx_fatal(const char *fmt, ...) __attribute__((noreturn, format(printf, 1, 2)));

void *sx_malloc(size_t n);
void *sx_calloc(size_t n, size_t sz);
void *sx_realloc(void *p, size_t n);

long sx_count_intdigits(long val, int sign);

/* Growable byte buffer used by all writers. */
typedef struct {
    char *data;
    size_t len, cap;
} sx_buf;

void sx_buf_init(sx_buf *b);
void sx_buf_free(sx_buf *b);
void sx_buf_reserve(sx_buf *b, size_t extra);
void sx_buf_putc(sx_buf *b, char c);
void sx_buf_puts(sx_buf *b, const char *s);
void sx_buf_printf(sx_buf *b, const char *fmt, ...) __attribute__((format(printf, 2, 3)));
void sx_buf_clear(sx_buf *b);
int sx_buf_write_file(const sx_buf *b, const char *path); /* creates/truncates */
int sx_buf_append_file(const sx_buf *b, FILE *f);

/* Verbosity-gated stdout logging. */
extern int sx_verbosity;
#define SX_LOG(level, ...) do { if (sx_verbosity >= (level)) { printf(__VA_ARGS__); } } while (0)

#endif
