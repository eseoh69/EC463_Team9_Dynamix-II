#include "md.h"
#include "cell_list.h"
#include "neighbor_list.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cell_mt.h"
#include <pthread.h>

//
// ======================================================
//  FULL O(N^2) FORCE KERNEL
// ======================================================
//
double md_compute_forces_full(Particle *p, const SimParams *sp)
{
    size_t N = sp->N;
    double eps = sp->epsilon;
    double sig = sp->sigma;
    double L   = sp->L;
    double rc2 = sp->rc2;

    const double tiny2 = 1e-30;
    

    long long npairs = 0;

    // Reset forces
    for (size_t i = 0; i < N; i++)
        p[i].fx = p[i].fy = p[i].fz = 0.0;

    double U = 0.0;

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = i + 1; j < N; j++)
        {
            double dx = p[i].x - p[j].x;
            double dy = p[i].y - p[j].y;
            double dz = p[i].z - p[j].z;

            md_minimage(&dx, L);
            md_minimage(&dy, L);
            md_minimage(&dz, L);

            
            double r2 = dx*dx + dy*dy + dz*dz;
            const double rmin2 = 0.2;  // rmin = 0.01
            if (r2 < rmin2) r2 = rmin2;
            if (r2 < tiny2 || r2 > rc2)
                continue;
            npairs++;

            

            double inv_r2 = 1.0 / r2;
            double sig2_r2 = (sig*sig) * inv_r2;
            double sr6  = sig2_r2 * sig2_r2 * sig2_r2;
            double sr12 = sr6 * sr6;

            double F_over_r = 24.0 * eps * (2.0*sr12 - sr6) * inv_r2;

            double fx = F_over_r * dx;
            double fy = F_over_r * dy;
            double fz = F_over_r * dz;

            p[i].fx += fx;
            p[i].fy += fy;
            p[i].fz += fz;

            p[j].fx -= fx;
            p[j].fy -= fy;
            p[j].fz -= fz;

            U += 4.0 * eps * (sr12 - sr6);
        }
    }

    //printf("min_r = %.6g\n", sqrt(min_r2));
    //printf("FULL npairs = %lld\n", npairs);
    return U;
}


//
// ======================================================
//  FULL O(N^2) VELOCITY VERLET INTEGRATOR
// ======================================================
//

static inline double wrap_pos(double x, double L)
    {
    x = fmod(x, L);
    if (x < 0) x += L;
    return x;
    }


double md_integrate_full(Particle *p, const SimParams *sp, double *Kout)
{
    size_t N = sp->N;
    double dt = sp->dt;
    double half = 0.5 * dt;
    double L = sp->L;

    // Half-step velocities
    for (size_t i = 0; i < N; i++) {
        p[i].vx += half * p[i].fx;
        p[i].vy += half * p[i].fy;
        p[i].vz += half * p[i].fz;
    }

    // Update positions
    for (size_t i = 0; i < N; i++) {
        p[i].x += dt * p[i].vx;
        p[i].y += dt * p[i].vy;
        p[i].z += dt * p[i].vz;

        // periodic wrap
        p[i].x = wrap_pos(p[i].x, L);
        p[i].y = wrap_pos(p[i].y, L);
        p[i].z = wrap_pos(p[i].z, L);
        if (!(p[i].x >= 0 && p[i].x < L &&
        p[i].y >= 0 && p[i].y < L &&
        p[i].z >= 0 && p[i].z < L)) {
        fprintf(stderr, "Position out of bounds after wrap: i=%zu x=%g y=%g z=%g L=%g\n",
            i, p[i].x, p[i].y, p[i].z, L);
        exit(1);
        }
    }

    // New forces
    double U = md_compute_forces_full(p, sp);

    // Final half-step + compute kinetic energy
    double K = 0.0;
    for (size_t i = 0; i < N; i++)
    {
        p[i].vx += half * p[i].fx;
        p[i].vy += half * p[i].fy;
        p[i].vz += half * p[i].fz;

        K += 0.5 * (p[i].vx*p[i].vx +
                    p[i].vy*p[i].vy +
                    p[i].vz*p[i].vz);
    }

    *Kout = K;
    return U;
}


/*

// -----------------------------
// Thread args for force compute
// -----------------------------
typedef struct {
    Particle *p;
    const SimParams *sp;
    size_t i_begin;
    size_t i_end;

    // Per-thread force accumulators (length N each)
    double *fx_local;
    double *fy_local;
    double *fz_local;

    // Per-thread potential energy
    double U_local;
} ForceThreadArgs;

static void *force_thread_fn(void *argp)
{

    

    ForceThreadArgs *a = (ForceThreadArgs *)argp;
    const size_t N = a->sp->N;

    // IMPORTANT: clear the per-thread force arrays every call
    memset(a->fx_local, 0, N * sizeof(double));
    memset(a->fy_local, 0, N * sizeof(double));
    memset(a->fz_local, 0, N * sizeof(double));

    Particle *p = a->p;
    const SimParams *sp = a->sp;

    //const size_t N = sp->N;
    const double eps = sp->epsilon;
    const double sig = sp->sigma;
    const double L   = sp->L;
    const double rc2 = sp->rc2;

    const double tiny2 = 1e-30;

    double U = 0.0;


    for (size_t i = a->i_begin; i < a->i_end; i++) {
        for (size_t j = i + 1; j < N; j++) {

            double dx = p[i].x - p[j].x;
            double dy = p[i].y - p[j].y;
            double dz = p[i].z - p[j].z;

            md_minimage(&dx, L);
            md_minimage(&dy, L);
            md_minimage(&dz, L);

            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < tiny2 || r2 > rc2) continue;

            double inv_r2   = 1.0 / r2;
            double sig2_r2  = (sig * sig) * inv_r2;
            double sr6      = sig2_r2 * sig2_r2 * sig2_r2;
            double sr12     = sr6 * sr6;

            double F_over_r = 24.0 * eps * (2.0*sr12 - sr6) * inv_r2;

            double fx = F_over_r * dx;
            double fy = F_over_r * dy;
            double fz = F_over_r * dz;

            // Accumulate into thread-private arrays to avoid races
            a->fx_local[i] += fx;
            a->fy_local[i] += fy;
            a->fz_local[i] += fz;

            a->fx_local[j] -= fx;
            a->fy_local[j] -= fy;
            a->fz_local[j] -= fz;

            U += 4.0 * eps * (sr12 - sr6);
        }
    }

    a->U_local = U;
    return NULL;
}


// ---------------------------------------------------------
// Pthreads force kernel (safe, deterministic up to FP order)
// ---------------------------------------------------------
double md_compute_forces_full_pthreads(Particle *p, const SimParams *sp, int nthreads)
{
    const size_t N = sp->N;

    if (nthreads < 1) nthreads = 1;
    if ((size_t)nthreads > N) nthreads = (int)N; // don’t spawn more threads than particles

    // Reset forces (serial; cheap relative to O(N^2))
    for (size_t i = 0; i < N; i++)
        p[i].fx = p[i].fy = p[i].fz = 0.0;

    // Allocate per-thread force buffers: nthreads * N doubles each component
    double *fx_all = (double *)calloc((size_t)nthreads * N, sizeof(double));
    double *fy_all = (double *)calloc((size_t)nthreads * N, sizeof(double));
    double *fz_all = (double *)calloc((size_t)nthreads * N, sizeof(double));
    if (!fx_all || !fy_all || !fz_all) {
        fprintf(stderr, "md_compute_forces_full_pthreads: allocation failed\n");
        free(fx_all); free(fy_all); free(fz_all);
        // Fallback to serial
        return md_compute_forces_full(p, sp);
    }

    pthread_t *threads = (pthread_t *)malloc((size_t)nthreads * sizeof(pthread_t));
    ForceThreadArgs *args = (ForceThreadArgs *)malloc((size_t)nthreads * sizeof(ForceThreadArgs));
    if (!threads || !args) {
        fprintf(stderr, "md_compute_forces_full_pthreads: allocation failed\n");
        free(fx_all); free(fy_all); free(fz_all);
        free(threads); free(args);
        return md_compute_forces_full(p, sp);
    }

    // Partition i range across threads
    for (int t = 0; t < nthreads; t++) {
        size_t i0 = (N * (size_t)t) / (size_t)nthreads;
        size_t i1 = (N * (size_t)(t + 1)) / (size_t)nthreads;

        args[t].p = p;
        args[t].sp = sp;
        args[t].i_begin = i0;
        args[t].i_end = i1;

        args[t].fx_local = fx_all + (size_t)t * N;
        args[t].fy_local = fy_all + (size_t)t * N;
        args[t].fz_local = fz_all + (size_t)t * N;

        args[t].U_local = 0.0;

        int rc = pthread_create(&threads[t], NULL, force_thread_fn, &args[t]);
        if (rc != 0) {
            // If thread creation fails, join what we started and fallback
            for (int k = 0; k < t; k++) pthread_join(threads[k], NULL);
            free(fx_all); free(fy_all); free(fz_all);
            free(threads); free(args);
            return md_compute_forces_full(p, sp);
        }
    }

    double U = 0.0;
    for (int t = 0; t < nthreads; t++) {
        pthread_join(threads[t], NULL);
        U += args[t].U_local;
    }

    // Reduce per-thread forces into p[i].f*
    for (size_t i = 0; i < N; i++) {
        double fx = 0.0, fy = 0.0, fz = 0.0;
        for (int t = 0; t < nthreads; t++) {
            size_t base = (size_t)t * N + i;
            fx += fx_all[base];
            fy += fy_all[base];
            fz += fz_all[base];
        }
        p[i].fx = fx;
        p[i].fy = fy;
        p[i].fz = fz;
    }

    free(fx_all); free(fy_all); free(fz_all);
    free(threads);
    free(args);

    return U;
}


// -----------------------------
// Thread args for integrate loops
// -----------------------------
typedef struct {
    Particle *p;
    const SimParams *sp;
    size_t begin;
    size_t end;
    double half_dt;
    double dt;
    double L;
    double K_local;
} IntegrateThreadArgs;

static void *integrate_halfstep_fn(void *argp)
{
    IntegrateThreadArgs *a = (IntegrateThreadArgs *)argp;
    Particle *p = a->p;

    for (size_t i = a->begin; i < a->end; i++) {
        p[i].vx += a->half_dt * p[i].fx;
        p[i].vy += a->half_dt * p[i].fy;
        p[i].vz += a->half_dt * p[i].fz;
    }
    return NULL;
}

static void *integrate_poswrap_fn(void *argp)
{
    IntegrateThreadArgs *a = (IntegrateThreadArgs *)argp;
    Particle *p = a->p;
    const double dt = a->dt;
    const double L  = a->L;

    for (size_t i = a->begin; i < a->end; i++) {
        p[i].x += dt * p[i].vx;
        p[i].y += dt * p[i].vy;
        p[i].z += dt * p[i].vz;

        // periodic wrap
        p[i].x = wrap_pos(p[i].x, L);
        p[i].y = wrap_pos(p[i].y, L);
        p[i].z = wrap_pos(p[i].z, L);
    
    }
    return NULL;
}

static void *integrate_finalhalf_K_fn(void *argp)
{
    IntegrateThreadArgs *a = (IntegrateThreadArgs *)argp;
    Particle *p = a->p;

    double K = 0.0;
    for (size_t i = a->begin; i < a->end; i++) {
        p[i].vx += a->half_dt * p[i].fx;
        p[i].vy += a->half_dt * p[i].fy;
        p[i].vz += a->half_dt * p[i].fz;

        K += 0.5 * (p[i].vx*p[i].vx +
                    p[i].vy*p[i].vy +
                    p[i].vz*p[i].vz);
    }
    a->K_local = K;
    return NULL;
}

static int spawn_partitioned_threads(
    pthread_t *threads,
    IntegrateThreadArgs *args,
    int nthreads,
    size_t N,
    void *(*fn)(void *),
    Particle *p,
    const SimParams *sp,
    double dt,
    double half_dt,
    double L
){
    for (int t = 0; t < nthreads; t++) {
        size_t b = (N * (size_t)t) / (size_t)nthreads;
        size_t e = (N * (size_t)(t + 1)) / (size_t)nthreads;

        args[t].p = p;
        args[t].sp = sp;
        args[t].begin = b;
        args[t].end = e;
        args[t].dt = dt;
        args[t].half_dt = half_dt;
        args[t].L = L;
        args[t].K_local = 0.0;

        int rc = pthread_create(&threads[t], NULL, fn, &args[t]);
        if (rc != 0) {
            for (int k = 0; k < t; k++) pthread_join(threads[k], NULL);
            return -1;
        }
    }
    for (int t = 0; t < nthreads; t++) pthread_join(threads[t], NULL);
    return 0;
}

// ---------------------------------------------------------
// Pthreads Velocity Verlet using pthread force kernel above
// ---------------------------------------------------------
double md_integrate_full_pthreads(Particle *p, const SimParams *sp, double *Kout, int nthreads)
{
    const size_t N = sp->N;
    if (nthreads < 1) nthreads = 1;
    if ((size_t)nthreads > N) nthreads = (int)N;

    const double dt   = sp->dt;
    const double half = 0.5 * dt;
    const double L    = sp->L;

    pthread_t *threads = (pthread_t *)malloc((size_t)nthreads * sizeof(pthread_t));
    IntegrateThreadArgs *args = (IntegrateThreadArgs *)malloc((size_t)nthreads * sizeof(IntegrateThreadArgs));
    if (!threads || !args) {
        free(threads); free(args);
        // fallback serial
        return md_integrate_full(p, sp, Kout);
    }

    // Half-step velocities (parallel)
    if (spawn_partitioned_threads(threads, args, nthreads, N,
                                  integrate_halfstep_fn, p, sp, dt, half, L) != 0) {
        free(threads); free(args);
        return md_integrate_full(p, sp, Kout);
    }

    // Update positions + periodic wrap (parallel)
    if (spawn_partitioned_threads(threads, args, nthreads, N,
                                  integrate_poswrap_fn, p, sp, dt, half, L) != 0) {
        free(threads); free(args);
        return md_integrate_full(p, sp, Kout);
    }

    // New forces (parallel O(N^2))
    double U = md_compute_forces_full_pthreads(p, sp, nthreads);

    // Final half-step + kinetic energy reduction (parallel)
    if (spawn_partitioned_threads(threads, args, nthreads, N,
                                  integrate_finalhalf_K_fn, p, sp, dt, half, L) != 0) {
        free(threads); free(args);
        return md_integrate_full(p, sp, Kout);
    }

    double K = 0.0;
    for (int t = 0; t < nthreads; t++) K += args[t].K_local;

    *Kout = K;

    free(threads);
    free(args);
    return U;
}
*/



// full_pthreads_owner.c
// O(N^2) "owner-computes" pthread implementation (race-free, simple, robust)


typedef struct {
    Particle *p;
    const SimParams *sp;
    size_t i0, i1;      // [i0, i1)
    double U_local;
} ForceArgs;

static void *force_owner_fn(void *argp)
{
    ForceArgs *a = (ForceArgs*)argp;
    Particle *p = a->p;
    const SimParams *sp = a->sp;

    const size_t N = sp->N;
    const double eps = sp->epsilon;
    const double sig = sp->sigma;
    const double L   = sp->L;
    const double rc2 = sp->rc2;
    const double tiny2 = 1e-30;

    double U = 0.0;

    for (size_t i = a->i0; i < a->i1; i++) {
        double fx = 0.0, fy = 0.0, fz = 0.0;

        const double xi = p[i].x, yi = p[i].y, zi = p[i].z;

        for (size_t j = 0; j < N; j++) {
            if (j == i) continue;

            double dx = xi - p[j].x;
            double dy = yi - p[j].y;
            double dz = zi - p[j].z;

            md_minimage(&dx, L);
            md_minimage(&dy, L);
            md_minimage(&dz, L);

            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < tiny2 || r2 > rc2) continue;

            double inv_r2  = 1.0 / r2;
            double sig2_r2 = (sig*sig) * inv_r2;
            double sr6     = sig2_r2 * sig2_r2 * sig2_r2;
            double sr12    = sr6 * sr6;

            double F_over_r = 24.0 * eps * (2.0*sr12 - sr6) * inv_r2;

            fx += F_over_r * dx;
            fy += F_over_r * dy;
            fz += F_over_r * dz;

            // Potential energy: each pair (i,j) is visited twice across all i,
            // so accumulate half to avoid double counting.
            U += 0.5 * 4.0 * eps * (sr12 - sr6);
        }

        p[i].fx = fx;
        p[i].fy = fy;
        p[i].fz = fz;
    }

    a->U_local = U;
    return NULL;
}

double md_compute_forces_full_pthreads(Particle *p, const SimParams *sp, int nthreads)
{
    const size_t N = sp->N;
    if (nthreads < 1) nthreads = 1;
    if ((size_t)nthreads > N) nthreads = (int)N;

    pthread_t *th = (pthread_t*)malloc((size_t)nthreads * sizeof(pthread_t));
    ForceArgs *a  = (ForceArgs*)malloc((size_t)nthreads * sizeof(ForceArgs));
    if (!th || !a) {
        free(th); free(a);
        return md_compute_forces_full(p, sp); // fallback
    }

    // partition i across threads
    for (int t = 0; t < nthreads; t++) {
        size_t i0 = (N * (size_t)t) / (size_t)nthreads;
        size_t i1 = (N * (size_t)(t + 1)) / (size_t)nthreads;

        a[t].p = p;
        a[t].sp = sp;
        a[t].i0 = i0;
        a[t].i1 = i1;
        a[t].U_local = 0.0;

        if (pthread_create(&th[t], NULL, force_owner_fn, &a[t]) != 0) {
            // join what started, fallback
            for (int k = 0; k < t; k++) pthread_join(th[k], NULL);
            free(th); free(a);
            return md_compute_forces_full(p, sp);
        }
    }

    double U = 0.0;
    for (int t = 0; t < nthreads; t++) {
        pthread_join(th[t], NULL);
        U += a[t].U_local;
    }

    free(th);
    free(a);
    return U;
}

// ------------------- Integrator (parallel phases) -------------------

typedef struct {
    Particle *p;
    const SimParams *sp;
    size_t i0, i1;
    double dt, half, L;
    double K_local;
} IntArgs;

static void *vv_halfstep_fn(void *argp)
{
    IntArgs *a = (IntArgs*)argp;
    Particle *p = a->p;
    const double half = a->half;

    for (size_t i = a->i0; i < a->i1; i++) {
        p[i].vx += half * p[i].fx;
        p[i].vy += half * p[i].fy;
        p[i].vz += half * p[i].fz;
    }
    return NULL;
}

static void *vv_poswrap_fn(void *argp)
{
    IntArgs *a = (IntArgs*)argp;
    Particle *p = a->p;
    const double dt = a->dt;
    const double L  = a->L;

    for (size_t i = a->i0; i < a->i1; i++) {
        p[i].x += dt * p[i].vx;
        p[i].y += dt * p[i].vy;
        p[i].z += dt * p[i].vz;

        p[i].x = wrap_pos(p[i].x, L);
        p[i].y = wrap_pos(p[i].y, L);
        p[i].z = wrap_pos(p[i].z, L);
    }
    return NULL;
}

static void *vv_finalhalf_K_fn(void *argp)
{
    IntArgs *a = (IntArgs*)argp;
    Particle *p = a->p;
    const double half = a->half;

    double K = 0.0;
    for (size_t i = a->i0; i < a->i1; i++) {
        p[i].vx += half * p[i].fx;
        p[i].vy += half * p[i].fy;
        p[i].vz += half * p[i].fz;

        K += 0.5 * (p[i].vx*p[i].vx + p[i].vy*p[i].vy + p[i].vz*p[i].vz);
    }
    a->K_local = K;
    return NULL;
}

static int run_partitioned(
    int nthreads, size_t N,
    void *(*fn)(void*),
    IntArgs *a, pthread_t *th,
    Particle *p, const SimParams *sp,
    double dt, double half, double L
){
    for (int t = 0; t < nthreads; t++) {
        size_t i0 = (N * (size_t)t) / (size_t)nthreads;
        size_t i1 = (N * (size_t)(t + 1)) / (size_t)nthreads;

        a[t].p = p;
        a[t].sp = sp;
        a[t].i0 = i0;
        a[t].i1 = i1;
        a[t].dt = dt;
        a[t].half = half;
        a[t].L = L;
        a[t].K_local = 0.0;

        if (pthread_create(&th[t], NULL, fn, &a[t]) != 0) {
            for (int k = 0; k < t; k++) pthread_join(th[k], NULL);
            return -1;
        }
    }
    for (int t = 0; t < nthreads; t++) pthread_join(th[t], NULL);
    return 0;
}

double md_integrate_full_pthreads(Particle *p, const SimParams *sp, double *Kout, int nthreads)
{
    const size_t N = sp->N;
    if (nthreads < 1) nthreads = 1;
    if ((size_t)nthreads > N) nthreads = (int)N;

    const double dt = sp->dt;
    const double half = 0.5 * dt;
    const double L = sp->L;

    pthread_t *th = (pthread_t*)malloc((size_t)nthreads * sizeof(pthread_t));
    IntArgs   *a  = (IntArgs*)malloc((size_t)nthreads * sizeof(IntArgs));
    if (!th || !a) {
        free(th); free(a);
        return md_integrate_full(p, sp, Kout);
    }

    if (run_partitioned(nthreads, N, vv_halfstep_fn, a, th, p, sp, dt, half, L) != 0) {
        free(th); free(a);
        return md_integrate_full(p, sp, Kout);
    }

    if (run_partitioned(nthreads, N, vv_poswrap_fn, a, th, p, sp, dt, half, L) != 0) {
        free(th); free(a);
        return md_integrate_full(p, sp, Kout);
    }

    double U = md_compute_forces_full_pthreads(p, sp, nthreads);

    if (run_partitioned(nthreads, N, vv_finalhalf_K_fn, a, th, p, sp, dt, half, L) != 0) {
        free(th); free(a);
        return md_integrate_full(p, sp, Kout);
    }

    double K = 0.0;
    for (int t = 0; t < nthreads; t++) K += a[t].K_local;

    *Kout = K;
    free(th);
    free(a);
    return U;
}





//
// ======================================================
//  CELL-LIST MD VELOCITY VERLET (wrapper)
//  Uses forces from md_compute_forces_cell()
// ======================================================
//


double md_integrate_cell(Particle *p, const SimParams *sp,
                         CellList *cl,
                         double *Kout)
{
    size_t N = sp->N;
    double dt = sp->dt;
    double half = 0.5 * dt;
    double L = sp->L;


    // Half-step velocities
    for (size_t i = 0; i < N; i++) {
        p[i].vx += half * p[i].fx;
        p[i].vy += half * p[i].fy;
        p[i].vz += half * p[i].fz;
    }

    // Update positions
    for (size_t i = 0; i < N; i++) {
        p[i].x += dt * p[i].vx;
        p[i].y += dt * p[i].vy;
        p[i].z += dt * p[i].vz;

        // periodic wrap
        p[i].x = wrap_pos(p[i].x, L);
        p[i].y = wrap_pos(p[i].y, L);
        p[i].z = wrap_pos(p[i].z, L);
        if (!(p[i].x >= 0 && p[i].x < L &&
        p[i].y >= 0 && p[i].y < L &&
        p[i].z >= 0 && p[i].z < L)) {
        fprintf(stderr, "Position out of bounds after wrap: i=%zu x=%g y=%g z=%g L=%g\n",
            i, p[i].x, p[i].y, p[i].z, L);
        exit(1);
}
    }

    // Must rebuild cell list every step (cheap)
    cell_list_build(cl, p, L);

    // New forces from cell-list LJ kernel
    double U = md_compute_forces_cell(p, sp, cl);

    // Final half-step + kinetic energy
    double K = 0.0;
    for (size_t i = 0; i < N; i++) {
        p[i].vx += half * p[i].fx;
        p[i].vy += half * p[i].fy;
        p[i].vz += half * p[i].fz;

        K += 0.5 * (p[i].vx*p[i].vx +
                    p[i].vy*p[i].vy +
                    p[i].vz*p[i].vz);
    }

    *Kout = K;
    return U;
}

// ======================================================
//  CELL-LIST MD VELOCITY VERLET (HALF-SHELL PARALLEL)
//  Uses forces from md_compute_forces_cell_mt
// ======================================================
double md_integrate_cell_mt(Particle *p, const SimParams *sp,
                                       CellList *cl,
                                       double *Kout)
{
    size_t N = sp->N;
    double dt = sp->dt;
    double half = 0.5 * dt;
    double L = sp->L;

    for (size_t i = 0; i < N; i++) {
        p[i].vx += half * p[i].fx;
        p[i].vy += half * p[i].fy;
        p[i].vz += half * p[i].fz;
    }

    for (size_t i = 0; i < N; i++) {
        p[i].x += dt * p[i].vx;
        p[i].y += dt * p[i].vy;
        p[i].z += dt * p[i].vz;

        p[i].x = wrap_pos(p[i].x, L);
        p[i].y = wrap_pos(p[i].y, L);
        p[i].z = wrap_pos(p[i].z, L);
    }

    cell_list_build(cl, p, L);

    // this function handles its own thread creation/joining internally
    double U = md_compute_forces_cell_mt(p, sp, cl);

    double K = 0.0;
    for (size_t i = 0; i < N; i++) {
        p[i].vx += half * p[i].fx;
        p[i].vy += half * p[i].fy;
        p[i].vz += half * p[i].fz;

        K += 0.5 * (p[i].vx*p[i].vx +
                    p[i].vy*p[i].vy +
                    p[i].vz*p[i].vz);
    }

    *Kout = K;
    return U;
}

