/*
#include "cell_list.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


//
// Helper: convert (ix,iy,iz) to flattened index
//
static inline int cell_index(int ix, int iy, int iz, int nc)
{
    return ix + nc * (iy + nc * iz);
}

//
// Helper: periodic wrap for integer cell coordinates
//
static inline int wrap(int a, int nc)
{
    if (a < 0)     return a + nc;
    if (a >= nc)  return a - nc;
    return a;
}

//
// Initialize cell list structure
// Store N inside the struct (per your choice 1B)
//
void cell_list_init(CellList *cl, size_t N, double L, double rc)
{
    cl->N = N;

    cl->nc = (int)(L / rc);
    if (cl->nc < 1) cl->nc = 1;

    cl->cell_size = L / cl->nc;
    cl->ncell     = cl->nc * cl->nc * cl->nc;

    // counts for each cell
    cl->counts = calloc(cl->ncell, sizeof(int));

    // generous capacity: up to 16 particles per cell (adjustable)
    cl->cells  = calloc(cl->ncell * 200, sizeof(int));
}

//
// Build the cell list for particle array p
//
void cell_list_build(CellList *cl, Particle *p, double L)
{
    int nc = cl->nc;
    size_t N = cl->N;

    // reset per-cell counts
    for (int c = 0; c < cl->ncell; c++)
        cl->counts[c] = 0;

    // clear cell storage
    size_t total_capacity = cl->N * 100;
    for (size_t idx = 0; idx < total_capacity; idx++)
        cl->cells[idx] = -1;

    double cs = cl->cell_size;

    // assign particles to cells
    for (size_t i = 0; i < N; i++)
    {
        int ix = (int)floor(p[i].x / cs);
        int iy = (int)floor(p[i].y / cs);
        int iz = (int)floor(p[i].z / cs);

        ix = wrap(ix, nc);
        iy = wrap(iy, nc);
        iz = wrap(iz, nc);

        int c = cell_index(ix, iy, iz, nc);

        int k = cl->counts[c]++;
        cl->cells[c * 200 + k] = (int)i;
    }
}

//
// Collect neighbors of particle i by scanning its 27 nearby cells.
// Used by BOTH: cell-list MD + neighbor list builder.
// Caller allocates buf[] of size at least N.
//
int cell_list_collect_neighbors(const CellList *cl,
                                size_t i,
                                const Particle *p,
                                double L,
                                int *buf)
{
    int nc = cl->nc;
    double cs = cl->cell_size;

    // Determine the home cell of particle i
    int ix = (int)floor(p[i].x / cs);
    int iy = (int)floor(p[i].y / cs);
    int iz = (int)floor(p[i].z / cs);

    ix = wrap(ix, nc);
    iy = wrap(iy, nc);
    iz = wrap(iz, nc);

    int count = 0;

    // Loop over the 27 surrounding cells (including home cell)
    for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dz = -1; dz <= 1; dz++)
    {
        int nx = wrap(ix + dx, nc);
        int ny = wrap(iy + dy, nc);
        int nz = wrap(iz + dz, nc);

        int c = cell_index(nx, ny, nz, nc);
        int base = c * 200;

        for (int k = 0; k < cl->counts[c]; k++) {
            int j = cl->cells[base + k];

            if (j < 0) continue;
            if ((size_t)j == i) continue;

            buf[count++] = j;
        }
    }

    return count;
}

// Collect indicies of neighbor cells for a given cell (including itself)
int cl_get_neighbor_cells(const CellList *cl, int cell_idx, int *neighbor_cells)
{
    int nc = cl->nc;

    // convert 1D index to 3D coordinates
    int z = cell_idx / (nc * nc);
    int rem = cell_idx % (nc * nc);
    int y = rem / nc;
    int x = rem % nc;

    int count = 0;

    // loop over the 3x3x3 neighborhood
    for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dz = -1; dz <= 1; dz++)
    {
        int nx = wrap(x + dx, nc);
        int ny = wrap(y + dy, nc);
        int nz = wrap(z + dz, nc);

        int neighbor_idx = cell_index(nx, ny, nz, nc);
        neighbor_cells[count++] = neighbor_idx;
    }

    return count;
}

//
// Full MD force computation using cell lists
// (per your choice 2C)
//
double md_compute_forces_cell(Particle *p, const SimParams *sp,
                              const CellList *cl)
{
    size_t N = sp->N;
    double eps = sp->epsilon;
    double sig = sp->sigma;
    double L   = sp->L;
    double rc2 = sp->rc2;

    // zero forces
    for (size_t i = 0; i < N; i++)
        p[i].fx = p[i].fy = p[i].fz = 0.0;

    double U = 0.0;
    const double tiny2 = 1e-30;

    // temp buffer to store neighbors
    int *buf = malloc(N * sizeof(int));

    for (size_t i = 0; i < N; i++)
    {
        int nnb = cell_list_collect_neighbors(cl, i, p, L, buf);

        for (int t = 0; t < nnb; t++) {
            int j = buf[t];
            if (j < 0) continue;

            double dx = p[i].x - p[j].x;
            double dy = p[i].y - p[j].y;
            double dz = p[i].z - p[j].z;

            md_minimage(&dx, L);
            md_minimage(&dy, L);
            md_minimage(&dz, L);

            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < tiny2 || r2 > rc2) continue;

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

            // Newton's third law
            p[j].fx -= fx;
            p[j].fy -= fy;
            p[j].fz -= fz;

            U += 4.0 * eps * (sr12 - sr6);
        }
    }

    free(buf);
    return U;
}

*/

#include "cell_list.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

//
// Helper: convert (ix,iy,iz) to flattened index
//
static inline int cell_index(int ix, int iy, int iz, int nc)
{
    return ix + nc * (iy + nc * iz);
}

//
// Helper: periodic wrap for integer cell coordinates
//
static inline int wrap(int a, int nc)
{
    a %= nc;
    if (a < 0) a += nc;
    return a;
}

// Grow the per-cell stride (row width) for cl->cells[]
static void cl_grow_stride(CellList *cl, int new_stride)
{
    int old_stride = cl->stride;
    if (new_stride <= old_stride) return;

    size_t new_total = (size_t)cl->ncell * (size_t)new_stride;
    int *new_cells = (int*)malloc(new_total * sizeof(int));
    if (!new_cells) {
        fprintf(stderr, "ERROR: cl_grow_stride malloc failed\n");
        exit(1);
    }

    // init to -1
    for (size_t idx = 0; idx < new_total; idx++) new_cells[idx] = -1;

    // copy old contents
    for (int c = 0; c < cl->ncell; c++) {
        int old_base = c * old_stride;
        int new_base = c * new_stride;
        int cnt = cl->counts[c];
        for (int k = 0; k < cnt; k++) {
            new_cells[new_base + k] = cl->cells[old_base + k];
        }
    }

    free(cl->cells);
    cl->cells = new_cells;
    cl->stride = new_stride;

    // now every cell can hold up to stride entries
    for (int c = 0; c < cl->ncell; c++) cl->cell_cap[c] = cl->stride;
}

//
// Initialize cell list structure
//
void cell_list_init(CellList *cl, size_t N, double L, double rc)
{
    cl->N = N;

    cl->nc = (int)(L / rc);
    if (cl->nc < 1) cl->nc = 1;

    cl->cell_size = L / cl->nc;
    while (cl->cell_size > rc && cl->nc < (int)L) {
    cl->nc++;
    cl->cell_size = L / cl->nc;
    }
    cl->ncell = cl->nc * cl->nc * cl->nc;

    // counts for each cell
    cl->counts = (int*)calloc(cl->ncell, sizeof(int));
    if (!cl->counts) { fprintf(stderr, "ERROR: calloc counts failed\n"); exit(1); }

    // start stride small; it will grow if needed
    cl->stride   = 16;
    cl->cell_cap = (int*)malloc((size_t)cl->ncell * sizeof(int));
    if (!cl->cell_cap) { fprintf(stderr, "ERROR: malloc cell_cap failed\n"); exit(1); }
    for (int c = 0; c < cl->ncell; c++) cl->cell_cap[c] = cl->stride;

    cl->cells = (int*)malloc((size_t)cl->ncell * (size_t)cl->stride * sizeof(int));
    if (!cl->cells) { fprintf(stderr, "ERROR: malloc cells failed\n"); exit(1); }

    // init to -1
    for (size_t idx = 0; idx < (size_t)cl->ncell * (size_t)cl->stride; idx++)
        cl->cells[idx] = -1;

    printf("nc=%d cell_size=%.6f rc=%.6f ratio=cell_size/rc=%.3f\n",
       cl->nc, cl->cell_size, rc, cl->cell_size/rc);
}

//
// Build the cell list for particle array p
//
void cell_list_build(CellList *cl, Particle *p, double L)
{
    int nc = cl->nc;
    size_t N = cl->N;

    // reset per-cell counts
    for (int c = 0; c < cl->ncell; c++)
        cl->counts[c] = 0;

    // clear cell storage (IMPORTANT: uses ncell*stride, not N*stride)
    size_t total_capacity = (size_t)cl->ncell * (size_t)cl->stride;
    for (size_t idx = 0; idx < total_capacity; idx++)
        cl->cells[idx] = -1;

    double cs = cl->cell_size;

    // assign particles to cells
    for (size_t i = 0; i < N; i++)
    {
        int ix = (int)floor(p[i].x / cs);
        int iy = (int)floor(p[i].y / cs);
        int iz = (int)floor(p[i].z / cs);

        ix = wrap(ix, nc);
        iy = wrap(iy, nc);
        iz = wrap(iz, nc);

        int c = cell_index(ix, iy, iz, nc);

        int k = cl->counts[c];

        // grow if this cell would overflow current stride
        if (k >= cl->cell_cap[c]) {
            int new_stride = cl->stride * 2;
            if (new_stride < k + 1) new_stride = k + 1;
            cl_grow_stride(cl, new_stride);

            // total_capacity changed; but we don't need to re-clear here because
            // we only grew due to insertion and cl_grow_stride initializes new slots to -1
        }

        cl->cells[c * cl->stride + k] = (int)i;
        cl->counts[c] = k + 1;
    }
}

//
// Collect neighbors of particle i by scanning its 27 nearby cells.
//
int cell_list_collect_neighbors(const CellList *cl,
                                size_t i,
                                const Particle *p,
                                double L,
                                int *buf)
{
    (void)L; // L not needed here; minimage handled in force routine
    int nc = cl->nc;
    double cs = cl->cell_size;

    // Determine the home cell of particle i
    int ix = (int)floor(p[i].x / cs);
    int iy = (int)floor(p[i].y / cs);
    int iz = (int)floor(p[i].z / cs);

    ix = wrap(ix, nc);
    iy = wrap(iy, nc);
    iz = wrap(iz, nc);

    int count = 0;

    // Loop over the 27 surrounding cells (including home cell)
    for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dz = -1; dz <= 1; dz++)
    {
        int nx = wrap(ix + dx, nc);
        int ny = wrap(iy + dy, nc);
        int nz = wrap(iz + dz, nc);

        int c = cell_index(nx, ny, nz, nc);
        int base = c * cl->stride;

        for (int k = 0; k < cl->counts[c]; k++) {
            int j = cl->cells[base + k];

            if (j < 0) continue;
            if ((size_t)j == i) continue;

            buf[count++] = j;
        }
    }

    return count;
}

// Collect indices of neighbor cells for a given cell (including itself)
int cl_get_neighbor_cells(const CellList *cl, int cell_idx, int *neighbor_cells)
{
    int nc = cl->nc;

    // convert 1D index to 3D coordinates
    int z = cell_idx / (nc * nc);
    int rem = cell_idx % (nc * nc);
    int y = rem / nc;
    int x = rem % nc;

    int count = 0;

    // loop over the 3x3x3 neighborhood
    for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dz = -1; dz <= 1; dz++)
    {
        int nx = wrap(x + dx, nc);
        int ny = wrap(y + dy, nc);
        int nz = wrap(z + dz, nc);

        int neighbor_idx = cell_index(nx, ny, nz, nc);
        neighbor_cells[count++] = neighbor_idx;
    }

    return count;
}

//
// Full MD force computation using cell lists
//
double md_compute_forces_cell(Particle *p, const SimParams *sp,
                              const CellList *cl)
{
    size_t N = sp->N;
    double eps = sp->epsilon;
    double sig = sp->sigma;
    double L   = sp->L;
    double rc2 = sp->rc2;

    // zero forces
    for (size_t i = 0; i < N; i++)
        p[i].fx = p[i].fy = p[i].fz = 0.0;

    double U = 0.0;
    const double tiny2 = 1e-30;

    long long npairs = 0;

    // temp buffer to store neighbors
    int *buf = (int*)malloc(N * sizeof(int));
    if (!buf) { fprintf(stderr, "ERROR: malloc buf failed\n"); exit(1); }

    for (size_t i = 0; i < N; i++)
    {
        int nnb = cell_list_collect_neighbors(cl, i, p, L, buf);

        for (int t = 0; t < nnb; t++) {
            int j = buf[t];
            if (j < 0) continue;

            // IMPORTANT: avoid double counting since we apply Newton's 3rd law
            if ((size_t)j <= i) continue;

            double dx = p[i].x - p[j].x;
            double dy = p[i].y - p[j].y;
            double dz = p[i].z - p[j].z;

            md_minimage(&dx, L);
            md_minimage(&dy, L);
            md_minimage(&dz, L);

            double r2 = dx*dx + dy*dy + dz*dz;
            const double rmin2 = 0.2;  // rmin = 0.01
            if (r2 < rmin2) r2 = rmin2;
            if (r2 < tiny2 || r2 > rc2) continue;

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

            // Newton's third law
            p[j].fx -= fx;
            p[j].fy -= fy;
            p[j].fz -= fz;

            U += 4.0 * eps * (sr12 - sr6);

        }
    }

    free(buf);
    //printf("CELL npairs = %lld\n", npairs);
    return U;
}
