#include "sx_tpool.h"
#include "sx_common.h"
#include <unistd.h>

typedef struct job { sx_job_fn fn; void *arg; struct job *next; } job;

struct sx_tpool {
    pthread_t *th; int n;
    pthread_mutex_t mu; pthread_cond_t cv_job, cv_done;
    job *head, *tail;
    int queued, running, stop;
};

int sx_hardware_threads(void)
{
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n > 0 ? (int)n : 1;
}

static job *pop_locked(sx_tpool *p)
{
    job *j = p->head;
    if (!j) return NULL;
    p->head = j->next; if (!p->head) p->tail = NULL;
    --p->queued; ++p->running;
    return j;
}

static void run_job(sx_tpool *p, job *j)
{
    j->fn(j->arg);
    free(j);
    pthread_mutex_lock(&p->mu);
    --p->running;
    if (p->running == 0 && p->queued == 0) pthread_cond_broadcast(&p->cv_done);
    pthread_mutex_unlock(&p->mu);
}

static void *worker(void *vp)
{
    sx_tpool *p = vp;
    for (;;) {
        pthread_mutex_lock(&p->mu);
        while (!p->head && !p->stop) pthread_cond_wait(&p->cv_job, &p->mu);
        if (!p->head && p->stop) { pthread_mutex_unlock(&p->mu); return NULL; }
        job *j = pop_locked(p);
        pthread_mutex_unlock(&p->mu);
        run_job(p, j);
    }
}

sx_tpool *sx_tpool_new(int nthreads)
{
    if (nthreads <= 0) nthreads = sx_hardware_threads();
    sx_tpool *p = sx_calloc(1, sizeof(*p));
    p->n = nthreads;
    pthread_mutex_init(&p->mu, NULL); pthread_cond_init(&p->cv_job, NULL); pthread_cond_init(&p->cv_done, NULL);
    p->th = sx_calloc((size_t)nthreads, sizeof(pthread_t));
    for (int i = 0; i < nthreads; ++i)
        if (pthread_create(&p->th[i], NULL, worker, p) != 0) sx_fatal("cannot create worker thread %d", i);
    return p;
}

int sx_tpool_size(const sx_tpool *p) { return p->n; }

void sx_tpool_submit(sx_tpool *p, sx_job_fn fn, void *arg)
{
    job *j = sx_malloc(sizeof(*j));
    j->fn = fn; j->arg = arg; j->next = NULL;
    pthread_mutex_lock(&p->mu);
    if (p->tail) p->tail->next = j; else p->head = j;
    p->tail = j; ++p->queued;
    pthread_cond_signal(&p->cv_job);
    pthread_mutex_unlock(&p->mu);
}

void sx_tpool_submit_front(sx_tpool *p, sx_job_fn fn, void *arg)
{
    job *j = sx_malloc(sizeof(*j));
    j->fn = fn; j->arg = arg;
    pthread_mutex_lock(&p->mu);
    j->next = p->head; p->head = j; if (!p->tail) p->tail = j;
    ++p->queued;
    pthread_cond_signal(&p->cv_job);
    pthread_mutex_unlock(&p->mu);
}

int sx_tpool_help(sx_tpool *p)
{
    pthread_mutex_lock(&p->mu);
    job *j = pop_locked(p);
    pthread_mutex_unlock(&p->mu);
    if (!j) return 0;
    run_job(p, j);
    return 1;
}

void sx_tpool_wait_all(sx_tpool *p)
{
    while (sx_tpool_help(p)) {}
    pthread_mutex_lock(&p->mu);
    while (p->queued > 0 || p->running > 0) pthread_cond_wait(&p->cv_done, &p->mu);
    pthread_mutex_unlock(&p->mu);
}

void sx_tpool_free(sx_tpool *p)
{
    if (!p) return;
    sx_tpool_wait_all(p);
    pthread_mutex_lock(&p->mu);
    p->stop = 1;
    pthread_cond_broadcast(&p->cv_job);
    pthread_mutex_unlock(&p->mu);
    for (int i = 0; i < p->n; ++i) pthread_join(p->th[i], NULL);
    pthread_mutex_destroy(&p->mu); pthread_cond_destroy(&p->cv_job); pthread_cond_destroy(&p->cv_done);
    free(p->th); free(p);
}
