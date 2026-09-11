#include "sx_common.h"

int MAX_IT = 10000;
int MAX_NAME = 100;
int MAX_LEAVES = 1000;
int NUM_BUFFER = 20;
int IO_BUFFER = 1000;
char TEST_CHAR = 7;
double FLOAT_PRECISION = 0.001;
int sx_verbosity = 1;

void sx_limits_from_env(void)
{
    const char *b;
    if ((b = getenv("SIMPHY_MAXIT")) != NULL) sscanf(b, "%u", (unsigned *)&MAX_IT);
    if ((b = getenv("SIMPHY_MAXNAME")) != NULL) sscanf(b, "%u", (unsigned *)&MAX_NAME);
    if ((b = getenv("SIMPHY_MAXLEAVES")) != NULL) sscanf(b, "%u", (unsigned *)&MAX_LEAVES);
    if ((b = getenv("SIMPHY_NUMBUFFER")) != NULL) sscanf(b, "%u", (unsigned *)&NUM_BUFFER);
    if ((b = getenv("SIMPHY_TESTCHAR")) != NULL) sscanf(b, "%c", &TEST_CHAR);
    if ((b = getenv("SIMPHY_IOBUFFER")) != NULL) sscanf(b, "%u", (unsigned *)&IO_BUFFER);
    if ((b = getenv("SIMPHY_FLOATPREC")) != NULL) sscanf(b, "%lf", &FLOAT_PRECISION);
}

void sx_check(long code, const char *msg)
{
    const char *m = msg ? msg : "";
    switch (code) {
    case SX_NO_ERROR:
        return;
    case SX_MEM_ERROR:
        fprintf(stderr, "\nDynamic memory error %s\n", m); break;
    case SX_IO_ERROR:
        fprintf(stderr, "\nIO_error %s\n", m); break;
    case SX_LOOP_ERROR:
        fprintf(stderr, "\nError: Infinite loop %s\n", m); break;
    case SX_SETTINGS_ERROR:
        fprintf(stderr, "\nSettings error: %s\nExecute simphy-x with the option -h in order to print the usage information\n", m); break;
    case SX_UNEXPECTED_VALUE:
        fprintf(stderr, "\nError: Unexpected value has been reached %s\n", m); break;
    case SX_DB_ERROR:
        fprintf(stderr, "\nError with the DB manipulation %s\n", m); break;
    case SX_TERMINATE_NOERROR:
        fflush(stdout); exit(EXIT_SUCCESS);
    default:
        return; /* pointer-like values: no error */
    }
    fflush(stdout);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

void sx_fatal(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fputc('\n', stderr);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
    fflush(stdout);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

void *sx_malloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (!p) sx_fatal("Dynamic memory error (malloc %zu bytes)", n);
    return p;
}
void *sx_calloc(size_t n, size_t sz)
{
    void *p = calloc(n ? n : 1, sz ? sz : 1);
    if (!p) sx_fatal("Dynamic memory error (calloc %zu x %zu bytes)", n, sz);
    return p;
}
void *sx_realloc(void *p, size_t n)
{
    void *q = realloc(p, n ? n : 1);
    if (!q) sx_fatal("Dynamic memory error (realloc %zu bytes)", n);
    return q;
}

long sx_count_intdigits(long val, int sign)
{
    long d = 1, c;
    if (val >= 0) for (c = 10; c <= val; c *= 10) d++;
    else for (c = -10; c >= val; c *= 10) d++;
    return (sign && c < 0) ? ++d : d;
}

/* ---- buffer ---- */
void sx_buf_init(sx_buf *b) { b->data = NULL; b->len = b->cap = 0; }
void sx_buf_free(sx_buf *b) { free(b->data); b->data = NULL; b->len = b->cap = 0; }
void sx_buf_clear(sx_buf *b) { b->len = 0; }
void sx_buf_reserve(sx_buf *b, size_t extra)
{
    if (b->len + extra + 1 > b->cap) {
        size_t nc = b->cap ? b->cap : 4096;
        while (nc < b->len + extra + 1) nc *= 2;
        b->data = sx_realloc(b->data, nc);
        b->cap = nc;
    }
}
void sx_buf_putc(sx_buf *b, char c)
{
    sx_buf_reserve(b, 1);
    b->data[b->len++] = c;
    b->data[b->len] = 0;
}
void sx_buf_puts(sx_buf *b, const char *s)
{
    size_t n = strlen(s);
    sx_buf_reserve(b, n);
    memcpy(b->data + b->len, s, n);
    b->len += n;
    b->data[b->len] = 0;
}
void sx_buf_printf(sx_buf *b, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    sx_buf_reserve(b, 64);
    int n = vsnprintf(b->data + b->len, b->cap - b->len, fmt, ap);
    va_end(ap);
    if (n < 0) sx_fatal("vsnprintf failure");
    if ((size_t)n >= b->cap - b->len) {
        sx_buf_reserve(b, (size_t)n + 1);
        va_start(ap, fmt);
        vsnprintf(b->data + b->len, b->cap - b->len, fmt, ap);
        va_end(ap);
    }
    b->len += (size_t)n;
}
int sx_buf_write_file(const sx_buf *b, const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) return SX_IO_ERROR;
    if (b->len && fwrite(b->data, 1, b->len, f) != b->len) { fclose(f); return SX_IO_ERROR; }
    fclose(f);
    return SX_NO_ERROR;
}
int sx_buf_append_file(const sx_buf *b, FILE *f)
{
    if (b->len && fwrite(b->data, 1, b->len, f) != b->len) return SX_IO_ERROR;
    return SX_NO_ERROR;
}
