#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>  // For clock_gettime

#include "config.h"
#include "md.h"
#include "md_io.h"
#include "pdb_importer.h"
#include "cell_list.h"
#include "neighbor_list.h"

#define MAX_ATOMS 50000

// ----------------------------------------------------
// Time measurement by clock_gettime()
// ----------------------------------------------------
double interval(struct timespec start, struct timespec end)
{
    struct timespec temp;
    temp.tv_sec = end.tv_sec - start.tv_sec;
    temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    if (temp.tv_nsec < 0) {
        temp.tv_sec = temp.tv_sec - 1;
        temp.tv_nsec = temp.tv_nsec + 1000000000;
    }
    return (((double)temp.tv_sec) + ((double)temp.tv_nsec) * 1.0e-9);
}

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

    printf("\n================ INITIALIZING SIMULATION ================\n");
    
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
           sp.rc, sp.dt, AUTOTUNE_N_TIMESTEPS);

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

    printf("\n================ BEGINNING AUTOTUNING ================");
    // ------------------------------------------------
    // 4. FULL O(N^2) MD - WITH TIMING
    // ------------------------------------------------
    printf("\n================ FULL O(N^2) MD ================\n");

    FILE *f_full = open_energy_csv("output/full_energies.csv");
    if (!f_full) {
        free(p0); free(p_full); free(p_cell); free(p_nbl);
        return 1;
    }

    struct timespec time_start, time_stop;
    
    // Start timing for FULL
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    // initial forces
    (void) md_compute_forces_full(p_full, &sp);

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_full(p_full, &sp, &K);
        fprintf(f_full, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    // End timing for FULL
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_full = interval(time_start, time_stop);

    fclose(f_full);
    io_write_pdb("output/full_positions.pdb", p_full, (size_t)N);
    printf("FULL MD done. Time: %.6f seconds\n", time_full);

    // ------------------------------------------------
    // 5. CELL-LIST MD - WITH TIMING
    // ------------------------------------------------
    printf("\n================ CELL-LIST MD ================\n");

    FILE *f_cell = open_energy_csv("output/cell_energies.csv");
    if (!f_cell) {
        free(p0); free(p_full); free(p_cell); free(p_nbl);
        return 1;
    }

    // Start timing for CELL-LIST
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    CellList cl;
    cell_list_init(&cl, (size_t)N, sp.L, sp.rc);
    cell_list_build(&cl, p_cell, sp.L);
    (void) md_compute_forces_cell(p_cell, &sp, &cl);

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_cell(p_cell, &sp, &cl, &K);
        fprintf(f_cell, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    // End timing for CELL-LIST
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_cell = interval(time_start, time_stop);

    fclose(f_cell);
    io_write_pdb("output/cell_positions.pdb", p_cell, (size_t)N);
    printf("CELL-LIST MD done. Time: %.6f seconds\n", time_cell);

    // ------------------------------------------------
    // 6. NEIGHBOR-LIST MD - WITH TIMING
    // ------------------------------------------------
    printf("\n================ NEIGHBOR-LIST MD ================\n");

    FILE *f_nbl = open_energy_csv("output/nbl_energies.csv");
    if (!f_nbl) {
        free(p0); free(p_full); free(p_cell); free(p_nbl);
        return 1;
    }

    // Start timing for NEIGHBOR-LIST
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

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

    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        double K, U;
        U = md_integrate_nbl(p_nbl, &sp, &nl, &cl2, &K);
        fprintf(f_nbl, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
    }

    // End timing for NEIGHBOR-LIST
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_nbl = interval(time_start, time_stop);

    fclose(f_nbl);
    io_write_pdb("output/nbl_positions.pdb", p_nbl, (size_t)N);
    printf("NEIGHBOR-LIST MD done. Time: %.6f seconds\n", time_nbl);

    // ------------------------------------------------
    // 7. COMPARE TIMES AND DETERMINE FASTEST
    // ------------------------------------------------
    printf("\n========================================\n");
    printf("PERFORMANCE COMPARISON\n");
    printf("========================================\n");
    printf("FULL O(N^2):     %.6f seconds (%.2fx speedup)\n", 
           time_full, time_full/time_full);
    printf("CELL-LIST:       %.6f seconds (%.2fx speedup)\n", 
           time_cell, time_full/time_cell);
    printf("NEIGHBOR-LIST:   %.6f seconds (%.2fx speedup)\n", 
           time_nbl, time_full/time_nbl);
    printf("========================================\n");

    // Determine fastest
    const char *fastest_name;
    int fastest_method;
    double fastest_time = time_full;
    
    fastest_method = 0;
    fastest_name = "FULL O(N^2)";
    
    if (time_cell < fastest_time) {
        fastest_time = time_cell;
        fastest_method = 1;
        fastest_name = "CELL-LIST";
    }
    if (time_nbl < fastest_time) {
        fastest_time = time_nbl;
        fastest_method = 2;
        fastest_name = "NEIGHBOR-LIST";
    }

    printf("\nFASTEST METHOD: %s (%.6f seconds)\n", fastest_name, fastest_time);

    // ------------------------------------------------
    // 8. RUN FULL SIMULATION WITH FASTEST METHOD
    // ------------------------------------------------
    printf("\n========================================\n");
    printf("RUNNING FULL SIMULATION WITH: %s\nfor %d TIMESTEPS\n", 
           fastest_name, USER_N_TIMESTEPS);
    printf("========================================\n");

    // Reset to initial conditions for final run
    Particle *p_final = malloc(N * sizeof(Particle));
    copy_particles(p_final, p0, (size_t)N);

    FILE *f_final = open_energy_csv("output/final_energies.csv");
    if (!f_final) {
        free(p0); free(p_full); free(p_cell); free(p_nbl); free(p_final);
        return 1;
    }

    // Start timing for final run
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    // ----------------------------------------------------
    // FIX: Reset INITIAL FORCES before FINAL SIMULATION
    // ----------------------------------------------------
    if (fastest_method == 0) {
        md_compute_forces_full(p_final, &sp);
    }
    else if (fastest_method == 1) {
        CellList tmp_cl;
        cell_list_init(&tmp_cl, (size_t)N, sp.L, sp.rc);
        cell_list_build(&tmp_cl, p_final, sp.L);
        md_compute_forces_cell(p_final, &sp, &tmp_cl);
        free(tmp_cl.counts);
        free(tmp_cl.cells);
    }
    else if (fastest_method == 2) {
        CellList tmp_cl;
        cell_list_init(&tmp_cl, (size_t)N, sp.L, sp.rc);
        cell_list_build(&tmp_cl, p_final, sp.L);

        NeighborList tmp_nl;
        nbl_init(&tmp_nl, (size_t)N, sp.rc, 0.3 * sp.rc);
        nbl_build(&tmp_nl, &tmp_cl, p_final, sp.L, sp.rc, (size_t)N);

        md_compute_forces_nbl(p_final, &sp, &tmp_nl);

        free(tmp_cl.counts);
        free(tmp_cl.cells);
        free(tmp_nl.nb);
        free(tmp_nl.nb_index);
        free(tmp_nl.prev);
    }

    // ----------------------------------------------------
    // FINAL RUN LOOP (4000 steps)
    // ----------------------------------------------------
    if (fastest_method == 0) {
        for (int step = 0; step < USER_N_TIMESTEPS; step++) {
            double K, U;
            U = md_integrate_full(p_final, &sp, &K);
            fprintf(f_final, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
        }
    }
    else if (fastest_method == 1) {
        CellList cl_final;
        cell_list_init(&cl_final, (size_t)N, sp.L, sp.rc);
        cell_list_build(&cl_final, p_final, sp.L);

        for (int step = 0; step < USER_N_TIMESTEPS; step++) {
            double K, U;
            U = md_integrate_cell(p_final, &sp, &cl_final, &K);
            fprintf(f_final, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
        }

        free(cl_final.counts);
        free(cl_final.cells);
    }
    else if (fastest_method == 2) {
        CellList cl_final;
        cell_list_init(&cl_final, (size_t)N, sp.L, sp.rc);
        cell_list_build(&cl_final, p_final, sp.L);

        NeighborList nl_final;
        nbl_init(&nl_final, (size_t)N, sp.rc, 0.3 * sp.rc);
        nbl_build(&nl_final, &cl_final, p_final, sp.L, sp.rc, (size_t)N);

        for (int step = 0; step < USER_N_TIMESTEPS; step++) {
            double K, U;
            U = md_integrate_nbl(p_final, &sp, &nl_final, &cl_final, &K);
            fprintf(f_final, "%d,%.10f,%.10f,%.10f\n", step, K, U, K+U);
        }

        free(cl_final.counts);
        free(cl_final.cells);
        free(nl_final.nb);
        free(nl_final.nb_index);
        free(nl_final.prev);
    }

    // End timing for final run
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    double time_final = interval(time_start, time_stop);

    fclose(f_final);
    io_write_pdb("output/final_positions.pdb", p_final, (size_t)N);

    printf("Final simulation complete. Time: %.6f seconds\n", time_final);
    printf("Output: output/final_energies.csv, output/final_positions.pdb\n");

    // ------------------------------------------------
    // Cleanup
    // ------------------------------------------------
    free(p0);
    free(p_full);
    free(p_cell);
    free(p_nbl);
    free(p_final);

    printf("\nAll simulations complete.\n");
    return 0;
}
