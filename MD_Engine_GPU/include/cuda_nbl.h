#pragma once

#include <cstdio>

#ifdef __cplusplus
extern "C" {
#endif

#include "md.h"

// Run neighbor-list MD on the GPU for nsteps steps.
// Writes step,K,U,E to f_energy (CSV) each step.
// Returns wall-clock elapsed time in seconds.
double cuda_nbl(const SimParams *sp, Particle *p, FILE *f_energy, int nsteps);

#ifdef __cplusplus
}
#endif
