#include "cell_list.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

static inline int cell_idx(int ix, int iy, int iz, int nc) {
    return ix + nc * (iy + nc * iz);
}

static inline int wrap(int a, int nc) {
    if (a < 0)   return a + nc;
    if (a >= nc) return a - nc;
    return a;
}

void cell_list_init(CellList *cl, size_t N, double L, double rc) {
    cl->N = N;
    cl->nc = (int)(L / rc);
    if (cl->nc < 1) cl->nc = 1;
    cl->cell_size = L / cl->nc;
    cl->ncell = cl->nc * cl->nc * cl->nc;
    cl->counts = (int *)calloc(cl->ncell, sizeof(int));
    cl->cells  = (int *)malloc((size_t)cl->ncell * CELL_CAP * sizeof(int));
}

void cell_list_build(CellList *cl, Particle *p, double L) {
    const int nc = cl->nc;
    const double cs = cl->cell_size;

    memset(cl->counts, 0, cl->ncell * sizeof(int));

    for (size_t i = 0; i < cl->N; i++) {
        int ix = wrap((int)floor(p[i].x / cs), nc);
        int iy = wrap((int)floor(p[i].y / cs), nc);
        int iz = wrap((int)floor(p[i].z / cs), nc);
        int c  = cell_idx(ix, iy, iz, nc);
        int k  = cl->counts[c];
        if (k < CELL_CAP) {
            cl->cells[c * CELL_CAP + k] = (int)i;
            cl->counts[c]++;
        }
    }
}

void cell_list_free(CellList *cl) {
    free(cl->counts);
    free(cl->cells);
    cl->counts = NULL;
    cl->cells  = NULL;
}
