#ifndef CONFIG_H
#define CONFIG_H

// -----------------------------
// Lennard-Jones (unitless)
// -----------------------------
#define CONF_SIGMA       1.0
#define CONF_EPSILON     1.0

// -----------------------------
// LJ Cutoff (in Angstroms)
// -----------------------------
#define CONF_CUTOFF      30

// -----------------------------
// Time step + steps
// -----------------------------
#define CONF_DELTA_T     0.005
#define AUTOTUNE_N_TIMESTEPS 2000
#define USER_N_TIMESTEPS     4000

// -----------------------------
// Box fallback
// -----------------------------
#define CONF_BOX_MIN     0.0
#define CONF_BOX_MAX     100.0

// -----------------------------
// GPU block size
// -----------------------------
#define GPU_BLOCK_SIZE   256

// -----------------------------
// Neighbor list parameters
// -----------------------------
#define NBL_SKIN_FACTOR  0.3     // skin = 0.3 * rc
#define MAX_NEIGHBORS    128     // max neighbors per particle

#endif
