#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "config.h"
#include "md.h"
#include "md_io.h"
#include "pdb_importer.h"
#include "cell_list.h"
#include "neighbor_list.h"

#define MAX_ATOMS 50000

// ----------------------------------------------------
// Helper: copy particle array
// ----------------------------------------------------
static void copy_particles(Particle *dst, const Particle *src, size_t N)
{
    for (size_t i = 0; i < N; i++) {
        dst[i] = src[i];
    }
}

// ----------------------------------------------------
// Helper: open energy CSV file
// ----------------------------------------------------
static FILE *open_energy_csv(const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) {
        perror("fopen energy csv");
        return NULL;
    }
    fprintf(f, "step,K,U,E\n");
    return f;
}

// ----------------------------------------------------
// main
// ----------------------------------------------------
int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("Usage: %s input.pdb\n", argv[0]);
        return 1;
    }

    const char *pdb_path = argv[1];

    // ------------------------------------------------
    // 1. Read PDB
    // ------------------------------------------------
    Particle *p0 = malloc(MAX_ATOMS * sizeof(Particle));
    if (!p0) {
        fprintf(stderr, "Allocation failure for p0\n");
        return 1;
    }

    double xmin, xmax, ymin, ymax, zmin, zmax;
    int N = pdb_importer(pdb_path, p0, MAX_ATOMS,
                         &xmin, &xmax,
                         &ymin, &ymax,
                         &zmin, &zmax);
    if (N <= 0) {
        fprintf(stderr, "Failed to read PDB: %s\n", pdb_path);
        free(p0);
        return 1;
    }

    printf("Loaded %d atoms from %s\n", N, pdb_path);
    printf("Bounds: dx = %.3f  dy = %.3f  dz = %.3f\n",
           xmax - xmin, ymax - ymin, zmax - zmin);

    // ------------------------------------------------
    // 2. Set up simulation parameters
    // ------------------------------------------------
    SimParams sp;
    sp.N       = (size_t)N;
    sp.sigma   = CONF_SIGMA;
    sp.epsilon = CONF_EPSILON;
    sp.rc      = CONF_CUTOFF;
    sp.rc2     = sp.rc * sp.rc;
    sp.dt      = CONF_DELTA_T;

    // Use largest span as box length
    sp.L = fmax(fmax(xmax - xmin, ymax - ymin), zmax - zmin);
    if (sp.L <= 0.0) sp.L = CONF_BOX_MAX;

    printf("Simulation box L = %.3f\n", sp.L);
    printf("rc = %.3f, dt = %.5f, steps = %d\n",
           sp.rc, sp.dt, CONF_N_TIMESTEPS);

    // ------------------------------------------------
    // 3. Make copies for each MD method
    // ------------------------------------------------
    Particle *p_full = malloc(N * sizeof(Particle));
    Particle *p_cell = malloc(N * sizeof(Particle));
    Particle *p_nbl  = malloc(N * sizeof(Particle));

    if (!p_full || !p_cell || !p_nbl) {
        fprintf(stderr, "Allocation failure for particle copies\n");
        free(p0);
        free(p_full); free(p_cell); free(p_nbl);
        return 1;
    }

    copy_particles(p_full, p0, (size_t)N);
    copy_particles(p_cell, p0, (size_t)N);
    copy_particles(p_nbl,  p0, (size_t)N);

    // ------------------------------------------------
    // 4. FULL O(N^2) MD
    // ------------------------------------------------
    printf("\n================ FULL O(N^2) MD ================\n");

    FILE *f_full = open_energy_csv("output/full_energies.csv");
    if (!f_full) {
        free(p0); free(p_full); free(p_cell); free(p_nbl);
        return 1;
    }

    // initial forces
    (void) md_compute_forces_full(p_full, &sp);

    for (int step = 0; step < CONF_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_full(p_full, &sp, &K);
        fprintf(f_full, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    fclose(f_full);
    io_write_pdb("output/full_positions.pdb", p_full, (size_t)N);
    printf("FULL MD done. Output: output/full_energies.csv, output/full_positions.pdb\n");

    // ------------------------------------------------
    // 5. CELL-LIST MD
    // ------------------------------------------------
    printf("\n================ CELL-LIST MD ================\n");

    FILE *f_cell = open_energy_csv("output/cell_energies.csv");
    if (!f_cell) {
        free(p0); free(p_full); free(p_cell); free(p_nbl);
        return 1;
    }

    CellList cl;
    cell_list_init(&cl, (size_t)N, sp.L, sp.rc);
    cell_list_build(&cl, p_cell, sp.L);
    (void) md_compute_forces_cell(p_cell, &sp, &cl);

    for (int step = 0; step < CONF_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_cell(p_cell, &sp, &cl, &K);
        fprintf(f_cell, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    fclose(f_cell);
    io_write_pdb("output/cell_positions.pdb", p_cell, (size_t)N);
    printf("CELL-LIST MD done. Output: output/cell_energies.csv, output/cell_positions.pdb\n");

    // ------------------------------------------------
    // 6. NEIGHBOR-LIST MD
    // ------------------------------------------------
    printf("\n================ NEIGHBOR-LIST MD ================\n");

    FILE *f_nbl = open_energy_csv("output/nbl_energies.csv");
    if (!f_nbl) {
        free(p0); free(p_full); free(p_cell); free(p_nbl);
        return 1;
    }

    // Cell list for NBL
    CellList cl2;
    cell_list_init(&cl2, (size_t)N, sp.L, sp.rc);
    cell_list_build(&cl2, p_nbl, sp.L);

    // Neighbor list
    NeighborList nl;
    nbl_init(&nl, (size_t)N, sp.rc, 0.3 * sp.rc);
    nbl_build(&nl, &cl2, p_nbl, sp.L, sp.rc, (size_t)N);

    // initial forces via NBL
    (void) md_compute_forces_nbl(p_nbl, &sp, &nl);

    for (int step = 0; step < CONF_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_nbl(p_nbl, &sp, &nl, &cl2, &K);
        fprintf(f_nbl, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    fclose(f_nbl);
    io_write_pdb("output/nbl_positions.pdb", p_nbl, (size_t)N);
    printf("NEIGHBOR-LIST MD done. Output: output/nbl_energies.csv, output/nbl_positions.pdb\n");

    // ------------------------------------------------
    // Cleanup
    // ------------------------------------------------
    free(p0);
    free(p_full);
    free(p_cell);
    free(p_nbl);

    printf("\nAll simulations complete.\n");
    return 0;
}
