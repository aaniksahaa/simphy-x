/* SPDX-License-Identifier: GPL-3.0-only */
/* Portability shims for Windows (see sx_port.h). Empty elsewhere. */
#include "sx_port.h"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#include <stdlib.h>

ssize_t sx_getline(char **lineptr, size_t *n, FILE *stream)
{
    if (lineptr == NULL || n == NULL || stream == NULL) return -1;
    if (*lineptr == NULL || *n == 0) {
        *n = 128;
        *lineptr = malloc(*n);
        if (*lineptr == NULL) return -1;
    }
    size_t len = 0;
    int c;
    while ((c = fgetc(stream)) != EOF) {
        if (len + 1 >= *n) {
            size_t grown = *n * 2;
            char *p = realloc(*lineptr, grown);
            if (p == NULL) return -1;
            *lineptr = p;
            *n = grown;
        }
        (*lineptr)[len++] = (char)c;
        if (c == '\n') break;
    }
    if (len == 0) return -1;
    (*lineptr)[len] = '\0';
    return (ssize_t)len;
}

void sx_port_init(void)
{
    _fmode = _O_BINARY;
    _setmode(_fileno(stdout), _O_BINARY);
    _setmode(_fileno(stderr), _O_BINARY);
}
#else
typedef int sx_port_unused; /* ISO C does not allow an empty file. */
#endif
