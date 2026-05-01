#pragma once

#include <stddef.h>
#include "cell_list.h"
#include "md.h"

// Manhattan (radius-1) cell-list force compute — owner-computes-i (serial).
// Each particle i accumulates its own force by visiting 7 Manhattan neighbor
// cells (self + ±x ±y ±z face neighbors). Pairs are visited twice so U is
// halved before returning. No threading, no half-shell.
//
// - p:  particle array (positions read, forces written)
// - sp: sim params (L, rc2, sigma, epsilon, N)
// - cl: cell list already built for current positions
// Returns potential energy U.
double md_compute_forces_manhat(
    Particle *p,
    const SimParams *sp,
    const CellList *cl
);

// Velocity-Verlet integrator using md_compute_forces_manhat.
// Returns potential energy U; writes kinetic energy to *K_out.
double md_integrate_manhat(
    Particle *p,
    const SimParams *sp,
    CellList *cl,
    double *K_out
);
