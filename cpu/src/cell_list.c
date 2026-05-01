#include "cell_list.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

static inline int cell_index(int ix, int iy, int iz, int nc)
{
    return ix + nc * (iy + nc * iz);
}

static inline int wrap(int a, int nc)
{
    a %= nc;
    if (a < 0) a += nc;
    return a;
}

// ─────────────────────────────────────────────────────────────────────────────
// cl_grow_stride
//
// Grows cl->cells from (ncell * old_stride) to (ncell * new_stride),
// copies existing contents, updates cl->stride and cl->cell_cap.
// ─────────────────────────────────────────────────────────────────────────────
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
    for (size_t idx = 0; idx < new_total; idx++)
        new_cells[idx] = -1;

    // copy old contents row by row
    for (int c = 0; c < cl->ncell; c++) {
        int cnt = cl->counts[c];
        for (int k = 0; k < cnt; k++)
            new_cells[(size_t)c * new_stride + k] =
                cl->cells[(size_t)c * old_stride + k];
    }

    free(cl->cells);
    cl->cells  = new_cells;
    cl->stride = new_stride;

    // every cell can now hold up to new_stride entries
    for (int c = 0; c < cl->ncell; c++)
        cl->cell_cap[c] = new_stride;
}

// ─────────────────────────────────────────────────────────────────────────────
// cell_list_init
// ─────────────────────────────────────────────────────────────────────────────
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

    // ── Pick an initial stride large enough to avoid growing on the first
    //    build.  Average occupancy = N/ncell; we give 4x headroom.
    //    Minimum of 16 to avoid tiny strides for small systems.
    int avg_occ = (cl->ncell > 0) ? (int)(N / (size_t)cl->ncell) + 1 : 1;
    cl->stride  = avg_occ * 4;
    if (cl->stride < 16) cl->stride = 16;

    cl->counts = (int*)calloc((size_t)cl->ncell, sizeof(int));
    if (!cl->counts) {
        fprintf(stderr, "ERROR: calloc counts failed\n"); exit(1);
    }

    cl->cell_cap = (int*)malloc((size_t)cl->ncell * sizeof(int));
    if (!cl->cell_cap) {
        fprintf(stderr, "ERROR: malloc cell_cap failed\n"); exit(1);
    }
    for (int c = 0; c < cl->ncell; c++)
        cl->cell_cap[c] = cl->stride;

    size_t total = (size_t)cl->ncell * (size_t)cl->stride;
    cl->cells = (int*)malloc(total * sizeof(int));
    if (!cl->cells) {
        fprintf(stderr, "ERROR: malloc cells failed\n"); exit(1);
    }
    for (size_t idx = 0; idx < total; idx++)
        cl->cells[idx] = -1;

    printf("nc=%d cell_size=%.6f rc=%.6f ratio=cell_size/rc=%.3f stride=%d\n",
           cl->nc, cl->cell_size, rc, cl->cell_size / rc, cl->stride);
}

// ─────────────────────────────────────────────────────────────────────────────
// cell_list_build
// ─────────────────────────────────────────────────────────────────────────────
void cell_list_build(CellList *cl, Particle *p, double L)
{
    int    nc = cl->nc;
    size_t N  = cl->N;

    // reset per-cell counts
    for (int c = 0; c < cl->ncell; c++)
        cl->counts[c] = 0;

    // clear cell storage using the CURRENT stride (may differ from init)
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

        // grow stride if this cell is full
        if (k >= cl->cell_cap[c]) {
            int new_stride = cl->stride * 2;
            if (new_stride < k + 1) new_stride = k + 1;
            cl_grow_stride(cl, new_stride);

            // After growing, total_capacity is stale — update it so the
            // NEXT call to cell_list_build clears the right amount.
            // (The current build continues correctly because cl_grow_stride
            //  already copied all previously written entries and the new
            //  slots are initialised to -1.)
            total_capacity = (size_t)cl->ncell * (size_t)cl->stride;
        }

        // safe to write: cl->stride is current after any growth above
        cl->cells[(size_t)c * (size_t)cl->stride + (size_t)k] = (int)i;
        cl->counts[c] = k + 1;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// cell_list_collect_neighbors
// ─────────────────────────────────────────────────────────────────────────────
int cell_list_collect_neighbors(const CellList *cl,
                                size_t i,
                                const Particle *p,
                                double L,
                                int *buf)
{
    (void)L;
    int    nc = cl->nc;
    double cs = cl->cell_size;

    int ix = (int)floor(p[i].x / cs);
    int iy = (int)floor(p[i].y / cs);
    int iz = (int)floor(p[i].z / cs);

    ix = wrap(ix, nc);
    iy = wrap(iy, nc);
    iz = wrap(iz, nc);

    int count = 0;

    for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dz = -1; dz <= 1; dz++)
    {
        int nx = wrap(ix + dx, nc);
        int ny = wrap(iy + dy, nc);
        int nz = wrap(iz + dz, nc);

        int c    = cell_index(nx, ny, nz, nc);
        int base = c * cl->stride;

        for (int k = 0; k < cl->counts[c]; k++) {
            int j = cl->cells[base + k];
            if (j < 0)           continue;
            if ((size_t)j == i)  continue;
            buf[count++] = j;
        }
    }

    return count;
}

// ─────────────────────────────────────────────────────────────────────────────
// cl_get_neighbor_cells
// ─────────────────────────────────────────────────────────────────────────────
int cl_get_neighbor_cells(const CellList *cl, int cell_idx, int *neighbor_cells)
{
    int nc = cl->nc;

    int z   = cell_idx / (nc * nc);
    int rem = cell_idx % (nc * nc);
    int y   = rem / nc;
    int x   = rem % nc;

    int count = 0;

    for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dz = -1; dz <= 1; dz++)
    {
        int nx = wrap(x + dx, nc);
        int ny = wrap(y + dy, nc);
        int nz = wrap(z + dz, nc);

        neighbor_cells[count++] = cell_index(nx, ny, nz, nc);
    }

    return count;
}

// ─────────────────────────────────────────────────────────────────────────────
// md_compute_forces_cell
// ─────────────────────────────────────────────────────────────────────────────
double md_compute_forces_cell(Particle *p, const SimParams *sp,
                              const CellList *cl)
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

    int *buf = (int*)malloc(N * sizeof(int));
    if (!buf) { fprintf(stderr, "ERROR: malloc buf failed\n"); exit(1); }

    for (size_t i = 0; i < N; i++)
    {
        int nnb = cell_list_collect_neighbors(cl, i, p, L, buf);

        for (int t = 0; t < nnb; t++) {
            int j = buf[t];
            if (j < 0)           continue;
            if ((size_t)j <= i)  continue;   // avoid double counting

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

    free(buf);
    return U;
}