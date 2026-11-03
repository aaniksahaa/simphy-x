/* SPDX-License-Identifier: GPL-3.0-only */
/* Minimal pthread thread pool with an unbounded FIFO queue. */
#ifndef SX_TPOOL_H
#define SX_TPOOL_H
#include <pthread.h>
typedef void (*sx_job_fn)(void *arg);
typedef struct sx_tpool sx_tpool;
sx_tpool *sx_tpool_new(int nthreads);           /* nthreads <= 0 → hardware threads */
int sx_tpool_size(const sx_tpool *p);
void sx_tpool_submit(sx_tpool *p, sx_job_fn fn, void *arg);
void sx_tpool_submit_front(sx_tpool *p, sx_job_fn fn, void *arg);   /* high priority */
/* Run one queued job on the calling thread if any; returns 1 if a job ran. */
int sx_tpool_help(sx_tpool *p);
void sx_tpool_wait_all(sx_tpool *p);            /* wait until queue empty and all jobs finished */
void sx_tpool_free(sx_tpool *p);
int sx_hardware_threads(void);
#endif
