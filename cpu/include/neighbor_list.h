#ifndef NEIGHBOR_LIST_H
#define NEIGHBOR_LIST_H

#include "md.h"
#include "cell_list.h"

//
// Proper named struct
//
typedef struct NeighborList {
    size_t N;
    int   *nb;
    int   *nb_index;
    Particle *prev;
    int total;
    int cap;
    double skin;
} NeighborList;

void nbl_init(NeighborList *nl, size_t N, double rc, double skin);
void nbl_copy_positions_N(NeighborList *nl, Particle *p, size_t N);
int  nbl_needs_rebuild(NeighborList *nl, Particle *p, size_t N, double L);

void nbl_build(NeighborList *nl,
               CellList *cl,
               Particle *p,
               double L,
               double rc,
               size_t N);

#endif
