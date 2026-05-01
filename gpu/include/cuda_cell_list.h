#pragma once

#include <cstdio>
#include <cuda_runtime.h>

// AoS particle layout shared by cell-list and neighbor-list (9 doubles)
struct Particle_cuda {
    double3 pos;
    double3 vel;
    double3 force;
};

#ifdef __cplusplus
extern "C" {
#endif

#include "md.h"

// Run cell-list MD on the GPU for nsteps steps.
// Writes step,K,U,E to f_energy (CSV) each step.
// Returns wall-clock elapsed time in seconds.
double cuda_cell_list(const SimParams *sp, Particle *p, FILE *f_energy, int nsteps);

#ifdef __cplusplus
}
#endif
