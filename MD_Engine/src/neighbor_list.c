#include "neighbor_list.h"
#include "cell_list.h"
#include "md.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

static inline double wrap_pos(double x, double L)
{
    x = fmod(x, L);
    if (x < 0) x += L;
    return x;
}

// ─────────────────────────────────────────────────────────────────────────────
// nbl_init
// ─────────────────────────────────────────────────────────────────────────────
void nbl_init(NeighborList *nl, size_t N, double rc, double skin)
{
    (void)rc;
    nl->N   = N;

    // Conservative starting capacity.
    // For large N with large skin the true count can be much larger than N*64,
    // so we start with N*128 and grow as needed.
    nl->cap = (int)(N * 128);
    if (nl->cap < 1024) nl->cap = 1024;

    nl->nb       = (int*)malloc((size_t)nl->cap * sizeof(int));
    nl->nb_index = (int*)malloc((N + 1) * sizeof(int));
    nl->prev     = (Particle*)malloc(N * sizeof(Particle));

    if (!nl->nb || !nl->nb_index || !nl->prev) {
        fprintf(stderr, "ERROR: nbl_init allocation failed\n");
        exit(1);
    }

    nl->total = 0;
    nl->skin  = skin;
}

// ─────────────────────────────────────────────────────────────────────────────
// nbl_ensure_capacity  (grows nl->nb if needed)
// ─────────────────────────────────────────────────────────────────────────────
static void nbl_ensure_capacity(NeighborList *nl, int needed_total)
{
    if (needed_total <= nl->cap) return;

    int newcap = nl->cap;
    while (newcap < needed_total) newcap *= 2;

    int *newnb = (int*)realloc(nl->nb, (size_t)newcap * sizeof(int));
    if (!newnb) {
        fprintf(stderr, "ERROR: realloc nb failed (needed %d)\n", newcap);
        exit(1);
    }
    nl->nb  = newnb;
    nl->cap = newcap;
}

// ─────────────────────────────────────────────────────────────────────────────
// Position snapshot helpers
// ─────────────────────────────────────────────────────────────────────────────
void nbl_copy_positions(NeighborList *nl, Particle *p)
{
    for (size_t i = 0; i < nl->N; i++) {
        nl->prev[i].x = p[i].x;
        nl->prev[i].y = p[i].y;
        nl->prev[i].z = p[i].z;
    }
}

void nbl_copy_positions_N(NeighborList *nl, Particle *p, size_t N)
{
    for (size_t i = 0; i < N; i++) {
        nl->prev[i].x = p[i].x;
        nl->prev[i].y = p[i].y;
        nl->prev[i].z = p[i].z;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// nbl_needs_rebuild
// ─────────────────────────────────────────────────────────────────────────────
int nbl_needs_rebuild(NeighborList *nl, Particle *p, size_t N, double L)
{
    double thresh  = 0.5 * nl->skin;
    double thresh2 = thresh * thresh;

    for (size_t i = 0; i < N; i++) {
        double dx = p[i].x - nl->prev[i].x;
        double dy = p[i].y - nl->prev[i].y;
        double dz = p[i].z - nl->prev[i].z;

        md_minimage(&dx, L);
        md_minimage(&dy, L);
        md_minimage(&dz, L);

        if (dx*dx + dy*dy + dz*dz > thresh2)
            return 1;
    }
    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
// nbl_build
//
// KEY FIX: call nbl_ensure_capacity(nl, nl->total + 1) before every write
// to nl->nb so we never overflow the buffer regardless of how many pairs
// are within (rc + skin).
// ─────────────────────────────────────────────────────────────────────────────
void nbl_build(NeighborList *nl, CellList *cl, Particle *p,
               double L, double rc, size_t N)
{
    double rc_skin  = rc + nl->skin;
    double rc_skin2 = rc_skin * rc_skin;

    nl->total = 0;
    for (size_t i = 0; i <= N; i++)
        nl->nb_index[i] = 0;

    int *buf = (int*)malloc(N * sizeof(int));
    if (!buf) {
        fprintf(stderr, "ERROR: nbl_build malloc buf failed\n");
        exit(1);
    }

    for (size_t i = 0; i < N; i++) {
        nl->nb_index[i] = nl->total;

        int nnb = cell_list_collect_neighbors(cl, i, p, L, buf);

        for (int t = 0; t < nnb; t++) {
            int j = buf[t];
            if (j < 0)           continue;
            if ((size_t)j <= i)  continue;   // half-shell: j > i only

            double dx = p[i].x - p[j].x;
            double dy = p[i].y - p[j].y;
            double dz = p[i].z - p[j].z;

            md_minimage(&dx, L);
            md_minimage(&dy, L);
            md_minimage(&dz, L);

            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 > rc_skin2) continue;

            // ── CAPACITY CHECK before every write ────────────────
            nbl_ensure_capacity(nl, nl->total + 1);
            nl->nb[nl->total++] = j;
        }
    }

    nl->nb_index[N] = nl->total;

    free(buf);
    nbl_copy_positions_N(nl, p, N);
}

// ─────────────────────────────────────────────────────────────────────────────
// md_compute_forces_nbl
// ─────────────────────────────────────────────────────────────────────────────
double md_compute_forces_nbl(Particle *p,
                             const SimParams *sp,
                             const NeighborList *nl)
{
    size_t N   = sp->N;
    double eps = sp->epsilon;
    double sig = sp->sigma;
    double L   = sp->L;
    double rc2 = sp->rc2;

    for (size_t i = 0; i < N; i++)
        p[i].fx = p[i].fy = p[i].fz = 0.0;

    double U = 0.0;
    const double tiny2 = 1e-30;
    const double rmin2 = 0.2;

    for (size_t i = 0; i < N; i++) {
        int start = nl->nb_index[i];
        int end   = nl->nb_index[i + 1];

        for (int k = start; k < end; k++) {
            int j = nl->nb[k];

            double dx = p[i].x - p[j].x;
            double dy = p[i].y - p[j].y;
            double dz = p[i].z - p[j].z;

            md_minimage(&dx, L);
            md_minimage(&dy, L);
            md_minimage(&dz, L);

            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < rmin2) r2 = rmin2;
            if (r2 < tiny2 || r2 > rc2) continue;

            double inv_r2  = 1.0 / r2;
            double sig2_r2 = (sig * sig) * inv_r2;
            double sr6     = sig2_r2 * sig2_r2 * sig2_r2;
            double sr12    = sr6 * sr6;

            double F_over_r = 24.0 * eps * (2.0*sr12 - sr6) * inv_r2;

            double fx = F_over_r * dx;
            double fy = F_over_r * dy;
            double fz = F_over_r * dz;

            p[i].fx += fx;  p[j].fx -= fx;
            p[i].fy += fy;  p[j].fy -= fy;
            p[i].fz += fz;  p[j].fz -= fz;

            U += 4.0 * eps * (sr12 - sr6);
        }
    }

    return U;
}

// ─────────────────────────────────────────────────────────────────────────────
// md_integrate_nbl
// ─────────────────────────────────────────────────────────────────────────────
double md_integrate_nbl(Particle *p,
                        const SimParams *sp,
                        NeighborList *nl,
                        CellList *cl,
                        double *Kout)
{
    size_t N    = sp->N;
    double dt   = sp->dt;
    double half = 0.5 * dt;
    double L    = sp->L;

    for (size_t i = 0; i < N; i++) {
        p[i].vx += half * p[i].fx;
        p[i].vy += half * p[i].fy;
        p[i].vz += half * p[i].fz;
    }

    for (size_t i = 0; i < N; i++) {
        p[i].x = wrap_pos(p[i].x + dt * p[i].vx, L);
        p[i].y = wrap_pos(p[i].y + dt * p[i].vy, L);
        p[i].z = wrap_pos(p[i].z + dt * p[i].vz, L);
    }

    if (nbl_needs_rebuild(nl, p, N, L)) {
        cell_list_build(cl, p, L);
        nbl_build(nl, cl, p, L, sp->rc, N);
    }

    double U = md_compute_forces_nbl(p, sp, nl);

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