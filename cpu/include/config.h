#ifndef CONFIG_H
#define CONFIG_H

// -----------------------------
// Lennard-Jones (reduced units)
// Positions must be in LJ units (sigma = 1). Use input/lj_1024.pdb.
// rho* = N*sigma^3/V ~ 0.8  ->  L ~ 10.86 for N=1024
// -----------------------------
#define CONF_SIGMA       1.0     // LJ sigma (reduced units)
#define CONF_EPSILON     1.0     // LJ epsilon (reduced units)

// -----------------------------
// LJ Cutoff: standard 2.5 sigma in LJ reduced units
// -----------------------------
#define CONF_CUTOFF      2.5

// -----------------------------
// Time step + steps
// -----------------------------
#define CONF_DELTA_T     0.005
#define AUTOTUNE_N_TIMESTEPS 1000

#define USER_N_TIMESTEPS 2000

// -----------------------------
// Box min/max (fallback — will be overwritten by PDB bounds)
// -----------------------------
#define CONF_BOX_MIN     0.0
#define CONF_BOX_MAX     11.0

#define N_THREADS        16

// -----------------------------
// Initial temperature (LJ units, kT = 1.0 is standard)
// -----------------------------
#define CONF_TEMPERATURE 1.0

#endif
