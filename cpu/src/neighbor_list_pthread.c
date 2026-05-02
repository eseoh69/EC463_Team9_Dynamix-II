/**
 * neighbor_list_pthread.c
 *
 * Half-Shell Neighbor List with Pthreads - LARGE SYSTEM OPTIMIZED
 * Thread count is now a runtime parameter passed into every public function.
 */

#include "neighbor_list_pthread.h"
#include "neighbor_list.h"
#include "cell_list.h"
#include "md.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

// ============================================================
//  MACOS BARRIER IMPLEMENTATION
//  (macOS doesn't have pthread_barrier_t)
// ============================================================

#ifdef __APPLE__

typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t  cond;
    int count;
    int tripCount;
} my_barrier_t;

static int my_barrier_init(my_barrier_t *barrier, unsigned int count)
{
    if (count == 0) return -1;
    if (pthread_mutex_init(&barrier->mutex, NULL) != 0) return -1;
    if (pthread_cond_init(&barrier->cond, NULL) != 0) {
        pthread_mutex_destroy(&barrier->mutex);
        return -1;
    }
    barrier->tripCount = count;
    barrier->count     = 0;
    return 0;
}

static int my_barrier_destroy(my_barrier_t *barrier)
{
    pthread_cond_destroy(&barrier->cond);
    pthread_mutex_destroy(&barrier->mutex);
    return 0;
}

static int my_barrier_wait(my_barrier_t *barrier)
{
    pthread_mutex_lock(&barrier->mutex);
    barrier->count++;
    if (barrier->count >= barrier->tripCount) {
        barrier->count = 0;
        pthread_cond_broadcast(&barrier->cond);
        pthread_mutex_unlock(&barrier->mutex);
        return 1;
    } else {
        pthread_cond_wait(&barrier->cond, &barrier->mutex);
        pthread_mutex_unlock(&barrier->mutex);
        return 0;
    }
}

#define barrier_t       my_barrier_t
#define barrier_init(b, count) my_barrier_init(b, count)
#define barrier_destroy(b)     my_barrier_destroy(b)
#define barrier_wait(b)        my_barrier_wait(b)

#else  // Linux / POSIX

#define barrier_t              pthread_barrier_t
#define barrier_init(b, count) pthread_barrier_init(b, NULL, count)
#define barrier_destroy(b)     pthread_barrier_destroy(b)
#define barrier_wait(b)        pthread_barrier_wait(b)

#endif

// ============================================================
//  HELPERS
// ============================================================

static inline double wrap_pos(double x, double L)
{
    x = fmod(x, L);
    if (x < 0) x += L;
    return x;
}

static void nbl_ensure_capacity(NeighborList *nl, int needed_total)
{
    if (needed_total <= nl->cap) return;
    int newcap = nl->cap > 0 ? nl->cap : 1024;
    while (newcap < needed_total) newcap *= 2;
    int *tmp = realloc(nl->nb, (size_t)newcap * sizeof(int));
    if (!tmp) {
        fprintf(stderr, "[NBL-PTHREAD] realloc nl->nb failed (need %d)\n", newcap);
        exit(EXIT_FAILURE);
    }
    nl->nb  = tmp;
    nl->cap = newcap;
}

// ============================================================
//  WORK MODES
// ============================================================

typedef enum {
    WORK_NONE = 0,
    WORK_COMPUTE_FORCES,
    WORK_REDUCE_FORCES,
    WORK_BUILD_NBL,
    WORK_SHUTDOWN
} WorkMode;

// ============================================================
//  THREAD DATA STRUCTURES
// ============================================================

typedef struct {
    int    tid;
    int    nthreads;
    size_t i_start;
    size_t i_end;

    Particle          *particles;
    const SimParams   *sp;
    const NeighborList *nl;
    const CellList    *cl;
    double  L;
    double  rc;
    double  rc_skin2;
    size_t  N;

    double *fx_priv;
    double *fy_priv;
    double *fz_priv;
    double  U_local;

    int *nb_local;
    int *nb_counts;
    int  nb_total;
    int  nb_cap;
} ThreadData;

// ============================================================
//  GLOBAL STATE  (fully dynamic — no fixed-size arrays)
// ============================================================

static struct {
    pthread_t   *threads;       /* [nthreads]   */
    ThreadData  *data;          /* [nthreads]   */

    barrier_t barrier_start;
    barrier_t barrier_end;

    volatile WorkMode work_mode;

    int    nthreads;            /* current pool size            */
    size_t alloc_N;             /* N the pool was allocated for */
    int    initialized;
    int    pool_running;
} g_state = {0};

// ============================================================
//  PERSISTENT WORKER THREAD
// ============================================================

static void *persistent_worker(void *arg)
{
    ThreadData *td = (ThreadData *)arg;

    while (1) {
        barrier_wait(&g_state.barrier_start);

        WorkMode mode = g_state.work_mode;

        if (mode == WORK_SHUTDOWN) break;

        /* ── FORCE COMPUTATION ───────────────────────────────── */
        if (mode == WORK_COMPUTE_FORCES) {
            const Particle     *p  = td->particles;
            const SimParams    *sp = td->sp;
            const NeighborList *nl = td->nl;

            const size_t N    = sp->N;
            const double eps  = sp->epsilon;
            const double sig  = sp->sigma;
            const double sig2 = sig * sig;
            const double L    = sp->L;
            const double rc2  = sp->rc2;
            const double tiny2 = 1e-30;

            double *fx = td->fx_priv;
            double *fy = td->fy_priv;
            double *fz = td->fz_priv;

            memset(fx, 0, N * sizeof(double));
            memset(fy, 0, N * sizeof(double));
            memset(fz, 0, N * sizeof(double));

            double U_local = 0.0;

            for (size_t i = td->i_start; i < td->i_end; i++) {
                const int nb_start = nl->nb_index[i];
                const int nb_end   = nl->nb_index[i + 1];

                const double xi = p[i].x;
                const double yi = p[i].y;
                const double zi = p[i].z;

                double fx_i = 0.0, fy_i = 0.0, fz_i = 0.0;

                for (int k = nb_start; k < nb_end; k++) {
                    const int j = nl->nb[k];

                    double dx = xi - p[j].x;
                    double dy = yi - p[j].y;
                    double dz = zi - p[j].z;

                    md_minimage(&dx, L);
                    md_minimage(&dy, L);
                    md_minimage(&dz, L);

                    double r2 = dx*dx + dy*dy + dz*dz;
                    if (r2 < tiny2 || r2 > rc2) continue;

                    double inv_r2   = 1.0 / r2;
                    double sig2_r2  = sig2 * inv_r2;
                    double sr6      = sig2_r2 * sig2_r2 * sig2_r2;
                    double sr12     = sr6 * sr6;
                    double F_over_r = 24.0 * eps * (2.0*sr12 - sr6) * inv_r2;

                    double fxij = F_over_r * dx;
                    double fyij = F_over_r * dy;
                    double fzij = F_over_r * dz;

                    fx_i  += fxij;  fy_i  += fyij;  fz_i  += fzij;
                    fx[j] -= fxij;  fy[j] -= fyij;  fz[j] -= fzij;

                    U_local += 4.0 * eps * (sr12 - sr6);
                }

                fx[i] += fx_i;
                fy[i] += fy_i;
                fz[i] += fz_i;
            }

            td->U_local = U_local;
        }

        /* ── PARALLEL REDUCTION ──────────────────────────────── */
        else if (mode == WORK_REDUCE_FORCES) {
            Particle *p        = td->particles;
            int       nthreads = td->nthreads;

            for (size_t i = td->i_start; i < td->i_end; i++) {
                double fx_sum = 0.0, fy_sum = 0.0, fz_sum = 0.0;
                for (int t = 0; t < nthreads; t++) {
                    fx_sum += g_state.data[t].fx_priv[i];
                    fy_sum += g_state.data[t].fy_priv[i];
                    fz_sum += g_state.data[t].fz_priv[i];
                }
                p[i].fx = fx_sum;
                p[i].fy = fy_sum;
                p[i].fz = fz_sum;
            }
        }

        /* ── NEIGHBOR LIST BUILD ─────────────────────────────── */
        else if (mode == WORK_BUILD_NBL) {
            const Particle *p        = td->particles;
            const CellList *cl       = td->cl;
            const double    L        = td->L;
            const double    rc_skin2 = td->rc_skin2;
            const size_t    N        = td->N;

            int *buf = malloc(N * sizeof(int));
            if (!buf) { td->nb_total = 0; goto build_done; }

            int local_idx = 0;

            for (size_t i = td->i_start; i < td->i_end; i++) {
                int count_i = 0;
                int ncand   = cell_list_collect_neighbors(cl, i, p, L, buf);

                for (int c = 0; c < ncand; c++) {
                    int j = buf[c];
                    if ((size_t)j <= i) continue;

                    double dx = p[i].x - p[j].x;
                    double dy = p[i].y - p[j].y;
                    double dz = p[i].z - p[j].z;

                    md_minimage(&dx, L);
                    md_minimage(&dy, L);
                    md_minimage(&dz, L);

                    double r2 = dx*dx + dy*dy + dz*dz;
                    if (r2 > rc_skin2) continue;

                    if (local_idx >= td->nb_cap) {
                        int  newcap = td->nb_cap * 2;
                        int *tmp    = realloc(td->nb_local,
                                              (size_t)newcap * sizeof(int));
                        if (!tmp) { td->nb_total = 0; free(buf); goto build_done; }
                        td->nb_local = tmp;
                        td->nb_cap   = newcap;
                    }
                    td->nb_local[local_idx++] = j;
                    count_i++;
                }

                td->nb_counts[i - td->i_start] = count_i;
            }

            td->nb_total = local_idx;
            free(buf);
            build_done:;
        }

        barrier_wait(&g_state.barrier_end);
    }

    return NULL;
}

// ============================================================
//  INITIALIZATION / CLEANUP
// ============================================================

void nbl_pthread_init(size_t N, int nthreads)
{
    /* Re-use existing pool if N and thread count are unchanged */
    if (g_state.initialized &&
        g_state.alloc_N  >= N &&
        g_state.nthreads == nthreads)
        return;

    if (g_state.initialized)
        nbl_pthread_cleanup();

    //printf("[NBL-PTHREAD] Initializing for N=%zu with %d threads\n",
           //N, nthreads);

    /* Allocate dynamic arrays */
    g_state.threads = malloc((size_t)nthreads * sizeof(pthread_t));
    g_state.data    = malloc((size_t)nthreads * sizeof(ThreadData));
    if (!g_state.threads || !g_state.data) {
        fprintf(stderr, "[NBL-PTHREAD] Failed to allocate thread arrays\n");
        exit(EXIT_FAILURE);
    }

    barrier_init(&g_state.barrier_start, nthreads + 1);
    barrier_init(&g_state.barrier_end,   nthreads + 1);

    size_t chunk = (N + (size_t)nthreads - 1) / (size_t)nthreads;

    for (int t = 0; t < nthreads; t++) {
        ThreadData *td = &g_state.data[t];

        td->tid      = t;
        td->nthreads = nthreads;
        td->i_start  = (size_t)t * chunk;
        td->i_end    = (size_t)(t + 1) * chunk;
        if (td->i_start > N) td->i_start = N;
        if (td->i_end   > N) td->i_end   = N;
        td->N = N;

        td->fx_priv = calloc(N, sizeof(double));
        td->fy_priv = calloc(N, sizeof(double));
        td->fz_priv = calloc(N, sizeof(double));
        if (!td->fx_priv || !td->fy_priv || !td->fz_priv) {
            fprintf(stderr, "[NBL-PTHREAD] Failed to allocate private arrays\n");
            exit(EXIT_FAILURE);
        }

        size_t range = (td->i_end > td->i_start)
                       ? (td->i_end - td->i_start) : 1;
        td->nb_cap   = (int)(range * 64);
        if (td->nb_cap < 1024) td->nb_cap = 1024;
        td->nb_local  = malloc((size_t)td->nb_cap * sizeof(int));
        td->nb_counts = malloc(range * sizeof(int));
        if (!td->nb_local || !td->nb_counts) {
            fprintf(stderr, "[NBL-PTHREAD] Failed to allocate build buffers\n");
            exit(EXIT_FAILURE);
        }
    }

    g_state.work_mode   = WORK_NONE;
    g_state.nthreads    = nthreads;
    g_state.alloc_N     = N;
    g_state.initialized = 1;

    for (int t = 0; t < nthreads; t++)
        pthread_create(&g_state.threads[t], NULL,
                       persistent_worker, &g_state.data[t]);
    g_state.pool_running = 1;

    double mem_mb = (3.0 * nthreads * N * sizeof(double)) / (1024.0 * 1024.0);
    //printf("[NBL-PTHREAD] Private array memory: %.2f MB\n", mem_mb);
    //printf("[NBL-PTHREAD] Persistent thread pool started\n");
}

void nbl_pthread_cleanup(void)
{
    if (!g_state.initialized) return;

    if (g_state.pool_running) {
        g_state.work_mode = WORK_SHUTDOWN;
        barrier_wait(&g_state.barrier_start);

        for (int t = 0; t < g_state.nthreads; t++)
            pthread_join(g_state.threads[t], NULL);
        g_state.pool_running = 0;
    }

    for (int t = 0; t < g_state.nthreads; t++) {
        ThreadData *td = &g_state.data[t];
        free(td->fx_priv);  free(td->fy_priv);  free(td->fz_priv);
        free(td->nb_local); free(td->nb_counts);
        td->fx_priv = td->fy_priv = td->fz_priv = NULL;
        td->nb_local = td->nb_counts = NULL;
    }

    free(g_state.threads); g_state.threads = NULL;
    free(g_state.data);    g_state.data    = NULL;

    barrier_destroy(&g_state.barrier_start);
    barrier_destroy(&g_state.barrier_end);

    g_state.initialized = 0;
    g_state.nthreads    = 0;
    g_state.alloc_N     = 0;

    printf("[NBL-PTHREAD] Cleanup complete\n");
}

// ============================================================
//  NEIGHBOR LIST BUILD
// ============================================================

void nbl_build_pthread(NeighborList *nl, CellList *cl, Particle *p,
                       double L, double rc, size_t N, int nthreads)
{
    nbl_pthread_init(N, nthreads);

    double rc_skin  = rc + nl->skin;
    double rc_skin2 = rc_skin * rc_skin;

    for (int t = 0; t < nthreads; t++) {
        g_state.data[t].particles = p;
        g_state.data[t].cl        = cl;
        g_state.data[t].L         = L;
        g_state.data[t].rc        = rc;
        g_state.data[t].rc_skin2  = rc_skin2;
        g_state.data[t].N         = N;
    }

    g_state.work_mode = WORK_BUILD_NBL;
    barrier_wait(&g_state.barrier_start);
    barrier_wait(&g_state.barrier_end);

    size_t chunk = (N + (size_t)nthreads - 1) / (size_t)nthreads;
    nl->total = 0;

    for (size_t i = 0; i < N; i++) {
        nl->nb_index[i] = nl->total;

        int owner = (int)(i / chunk);
        if (owner >= nthreads) owner = nthreads - 1;

        ThreadData *td = &g_state.data[owner];
        if (i >= td->i_start && i < td->i_end)
            nl->total += td->nb_counts[i - td->i_start];
    }
    nl->nb_index[N] = nl->total;

    nbl_ensure_capacity(nl, nl->total);

    for (int t = 0; t < nthreads; t++) {
        ThreadData *td  = &g_state.data[t];
        int src_offset  = 0;
        for (size_t i = td->i_start; i < td->i_end; i++) {
            int count = td->nb_counts[i - td->i_start];
            int dst   = nl->nb_index[i];
            memcpy(&nl->nb[dst], &td->nb_local[src_offset],
                   count * sizeof(int));
            src_offset += count;
        }
    }

    nbl_copy_positions_N(nl, p, N);
}

// ============================================================
//  FORCE COMPUTATION
// ============================================================

double md_compute_forces_nbl_pthread(Particle *p,
                                     const SimParams *sp,
                                     const NeighborList *nl,
                                     int nthreads)
{
    const size_t N = sp->N;

    nbl_pthread_init(N, nthreads);

    for (int t = 0; t < nthreads; t++) {
        g_state.data[t].particles = p;
        g_state.data[t].sp        = sp;
        g_state.data[t].nl        = nl;
        g_state.data[t].N         = N;
    }

    g_state.work_mode = WORK_COMPUTE_FORCES;
    barrier_wait(&g_state.barrier_start);
    barrier_wait(&g_state.barrier_end);

    g_state.work_mode = WORK_REDUCE_FORCES;
    barrier_wait(&g_state.barrier_start);
    barrier_wait(&g_state.barrier_end);

    double U_total = 0.0;
    for (int t = 0; t < nthreads; t++)
        U_total += g_state.data[t].U_local;

    return U_total;
}

// ============================================================
//  INTEGRATOR
// ============================================================

double md_integrate_nbl_pthread(Particle *p,
                                const SimParams *sp,
                                NeighborList *nl,
                                CellList *cl,
                                double *Kout,
                                int nthreads)
{
    const size_t N      = sp->N;
    const double dt     = sp->dt;
    const double half_dt = 0.5 * dt;
    const double L      = sp->L;

    for (size_t i = 0; i < N; i++) {
        p[i].vx += half_dt * p[i].fx;
        p[i].vy += half_dt * p[i].fy;
        p[i].vz += half_dt * p[i].fz;
    }

    for (size_t i = 0; i < N; i++) {
        p[i].x = wrap_pos(p[i].x + dt * p[i].vx, L);
        p[i].y = wrap_pos(p[i].y + dt * p[i].vy, L);
        p[i].z = wrap_pos(p[i].z + dt * p[i].vz, L);
    }

    if (nbl_needs_rebuild(nl, p, N, L)) {
        cell_list_build(cl, p, L);
        nbl_build_pthread(nl, cl, p, L, sp->rc, N, nthreads);
    }

    double U = md_compute_forces_nbl_pthread(p, sp, nl, nthreads);

    double K = 0.0;
    for (size_t i = 0; i < N; i++) {
        p[i].vx += half_dt * p[i].fx;
        p[i].vy += half_dt * p[i].fy;
        p[i].vz += half_dt * p[i].fz;
        K += 0.5 * (p[i].vx*p[i].vx + p[i].vy*p[i].vy + p[i].vz*p[i].vz);
    }

    *Kout = K;
    return U;
}

// ============================================================
//  UTILITIES
// ============================================================

int nbl_pthread_get_num_threads(void)
{
    return g_state.nthreads;
}


void nbl_pthread_print_info(void)
{
    printf("========================================================\n");
    printf("  NEIGHBOR-LIST PTHREAD (Half-Shell) - LARGE OPTIMIZED\n");
    printf("========================================================\n");
    printf("  Threads:              %d (runtime)\n", g_state.nthreads);
    printf("  Method:               Half-shell (j > i)\n");
    printf("  Race prevention:      Private force arrays\n");
    printf("  Reduction:            PARALLEL\n");
    printf("  Thread pool:          PERSISTENT (barrier sync)\n");
#ifdef __APPLE__
    printf("  Platform:             macOS (custom barrier)\n");
#else
    printf("  Platform:             Linux/POSIX\n");
#endif
    printf("========================================================\n");
}

size_t nbl_pthread_memory_usage(size_t N, int nthreads)
{
    return 3 * N * sizeof(double) * (size_t)nthreads;
}