#pragma once

#include <cstdio>
#include "md.h"
#include "cell_list.h"

#ifdef __cplusplus
extern "C" {
#endif

// Low-level GPU interface (CPU must build the CellList first)
int    gpu_manhattan_cuda_init(void);
void   gpu_manhattan_cuda_cleanup(void);
double gpu_manhattan_cuda_forces_and_U(Particle *p, const SimParams *sp, const CellList *cl);

// High-level interface matching cuda_full / cuda_cell_list / cuda_nbl
double cuda_manhattan(const SimParams *sp, Particle *p, FILE *f_energy, int nsteps);

#ifdef __cplusplus
}
#endif
