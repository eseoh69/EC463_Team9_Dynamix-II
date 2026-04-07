#ifndef MD_H
#define MD_H

#include <stddef.h>

//
// Particle structure
//
typedef struct {
    double x, y, z;
    double vx, vy, vz;
    double fx, fy, fz;
} Particle;

//
// Simulation parameters
//
typedef struct {
    size_t N;
    double L;

    double sigma;
    double epsilon;

    double rc;
    double rc2;

    double dt;
} SimParams;

//
// Minimum-image periodic BC
//
static inline void md_minimage(double *d, double L)
{
    if (*d >  0.5*L) *d -= L;
    if (*d < -0.5*L) *d += L;
}

#endif
