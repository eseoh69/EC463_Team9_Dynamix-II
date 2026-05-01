#include "manhat.h"

#include <math.h>
#include <stdlib.h>

#ifndef CELL_CAP
#define CELL_CAP 16
#endif

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static inline int wrap(int a, int nc) {
    if (a < 0) return a + nc;
    if (a >= nc) return a - nc;
    return a;
}

static inline int cell_index(int ix, int iy, int iz, int nc) {
    return ix + nc * (iy + nc * iz);
}

// Build the full 27-cell (3x3x3) neighborhood for the cell containing (x,y,z).
// 7-cell face-only was incorrect: edge/corner-adjacent cells share an edge/vertex
// and can have particles within rc that were missed, causing energy explosion.
static void neighbors27_for_pos(
    double x, double y, double z,
    double L, int nc,
    int out27[27]
) {
    double cs = L / (double)nc;

    int cx = (int)floor(x / cs);
    int cy = (int)floor(y / cs);
    int cz = (int)floor(z / cs);

    if (cx >= nc) cx = nc - 1;
    if (cy >= nc) cy = nc - 1;
    if (cz >= nc) cz = nc - 1;
    if (cx < 0)   cx = 0;
    if (cy < 0)   cy = 0;
    if (cz < 0)   cz = 0;

    int k = 0;
    for (int dz = -1; dz <= 1; dz++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dx = -1; dx <= 1; dx++)
        out27[k++] = cell_index(wrap(cx+dx,nc), wrap(cy+dy,nc), wrap(cz+dz,nc), nc);
}

static inline double wrap_pos(double x, double L) {
    x = fmod(x, L);
    if (x < 0) x += L;
    return x;
}

// ---------------------------------------------------------------------------
// Force compute — owner-computes-i, serial
// ---------------------------------------------------------------------------

double md_compute_forces_manhat(
    Particle *p,
    const SimParams *sp,
    const CellList *cl
) {
    const size_t N   = sp->N;
    const double L   = sp->L;
    const double rc2 = sp->rc2;
    const double sig = sp->sigma;
    const double eps = sp->epsilon;
    const int    nc  = cl->nc;

    // Zero forces
    for (size_t i = 0; i < N; i++) {
        p[i].fx = 0.0;
        p[i].fy = 0.0;
        p[i].fz = 0.0;
    }

    double U = 0.0;

    for (size_t i = 0; i < N; i++) {
        const double xi = p[i].x;
        const double yi = p[i].y;
        const double zi = p[i].z;

        int neigh[27];
        neighbors27_for_pos(xi, yi, zi, L, nc, neigh);

        double fxi = 0.0, fyi = 0.0, fzi = 0.0;
        double Ui  = 0.0;

        for (int n = 0; n < 27; n++) {
            const int c    = neigh[n];
            const int cnt  = cl->counts[c];
            const int base = c * CELL_CAP;

            for (int k = 0; k < cnt; k++) {
                const int j = cl->cells[base + k];
                if (j < 0 || (size_t)j == i) continue;

                double dx = xi - p[j].x;
                double dy = yi - p[j].y;
                double dz = zi - p[j].z;

                // minimum image convention
                dx -= L * nearbyint(dx / L);
                dy -= L * nearbyint(dy / L);
                dz -= L * nearbyint(dz / L);

                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 > rc2 || r2 < 1e-20) continue;

                const double inv_r2  = 1.0 / r2;
                const double sig2_r2 = (sig * sig) * inv_r2;
                const double sr6     = sig2_r2 * sig2_r2 * sig2_r2;
                const double sr12    = sr6 * sr6;

                Ui += 4.0 * eps * (sr12 - sr6);

                const double F_over_r = 24.0 * eps * (2.0 * sr12 - sr6) * inv_r2;
                fxi += F_over_r * dx;
                fyi += F_over_r * dy;
                fzi += F_over_r * dz;
            }
        }

        p[i].fx = fxi;
        p[i].fy = fyi;
        p[i].fz = fzi;
        U += Ui;
    }

    // Each pair (i,j) was visited twice: once from i's loop, once from j's.
    U *= 0.5;
    return U;
}

// ---------------------------------------------------------------------------
// Velocity-Verlet integrator
// ---------------------------------------------------------------------------

double md_integrate_manhat(
    Particle *p,
    const SimParams *sp,
    CellList *cl,
    double *K_out
) {
    const size_t N  = sp->N;
    const double dt = sp->dt;
    const double L  = sp->L;

    // Half-step velocity + full position update
    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;

        p[i].x += dt * p[i].vx;
        p[i].y += dt * p[i].vy;
        p[i].z += dt * p[i].vz;

        p[i].x = wrap_pos(p[i].x, L);
        p[i].y = wrap_pos(p[i].y, L);
        p[i].z = wrap_pos(p[i].z, L);
    }

    // Rebuild cell list for new positions
    cell_list_build(cl, p, L);

    // New forces
    double U = md_compute_forces_manhat(p, sp, cl);

    // Second half-step velocity update
    double K = 0.0;
    for (size_t i = 0; i < N; i++) {
        p[i].vx += 0.5 * dt * p[i].fx;
        p[i].vy += 0.5 * dt * p[i].fy;
        p[i].vz += 0.5 * dt * p[i].fz;

        K += 0.5 * (p[i].vx*p[i].vx + p[i].vy*p[i].vy + p[i].vz*p[i].vz);
    }

    if (K_out) *K_out = K;
    return U;
}
