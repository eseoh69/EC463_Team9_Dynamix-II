#ifndef CELL_LIST_H
#define CELL_LIST_H

#include "md.h"

#ifdef __cplusplus
extern "C" {
#endif

// CELL_CAP must match the value used in gpu_manhattan_cell_cap.cu
#define CELL_CAP 16

typedef struct CellList {
    size_t N;
    int nc;
    int ncell;
    double cell_size;
    int *counts;
    int *cells;
} CellList;

void cell_list_init(CellList *cl, size_t N, double L, double rc);
void cell_list_build(CellList *cl, Particle *p, double L);
void cell_list_free(CellList *cl);

#ifdef __cplusplus
}
#endif

#endif
