/* SPDX-License-Identifier: GPL-3.0-only */
/*
 * Portability shims for building SimPhy-X with MinGW-w64 on Windows. On other
 * systems this header is empty. Include it after the system headers.
 */
#ifndef SX_PORT_H
#define SX_PORT_H

#ifdef _WIN32
#include <stdio.h>
#include <sys/types.h>
#include <direct.h>

/* The Windows mkdir takes no permission bits. */
#define mkdir(path, mode) _mkdir(path)

/* The Windows C runtime has no getline. */
ssize_t sx_getline(char **lineptr, size_t *n, FILE *stream);
#define getline sx_getline

/* Opens every file in binary mode, so the output files have the same bytes
 * (no CRLF line ends) as on Linux. Called first in main. */
void sx_port_init(void);
#else
static inline void sx_port_init(void) {}
#endif

#endif
