#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "config.h"
#include "md.h"
#include "md_io.h"
#include "cuda_cell_list.h"
#include "cuda_nbl.h"
#include "cuda_full.h"
#include "gpu_manhattan_cell_cap.h"

static FILE *open_energy_csv(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) { perror("fopen energy csv"); return NULL; }
    fprintf(f, "step,K,U,E\n");
    return f;
}

static Particle *copy_particles(const Particle *src, size_t N) {
    Particle *dst = (Particle *)malloc(N * sizeof(Particle));
    if (dst) memcpy(dst, src, N * sizeof(Particle));
    return dst;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage: %s input_path\n", argv[0]);
        return 1;
    }

    const char *input_path = argv[1];

    printf("\n================ INITIALIZING SIMULATION ================\n");

    Particle *p0;
    SimParams sp;
    if (read_input(input_path, &p0, &sp) != 0) return 1;

    sp.sigma   = CONF_SIGMA;
    sp.epsilon = CONF_EPSILON;
    sp.dt      = CONF_DELTA_T;
    sp.rc      = fmin(0.5 * sp.L, CONF_CUTOFF);
    sp.rc2     = sp.rc * sp.rc;

    printf("Simulation box L = %.3f\n", sp.L);
    printf("rc = %.3f, dt = %.5f, steps = %d\n", sp.rc, sp.dt, AUTOTUNE_N_TIMESTEPS);

    Particle *p_full = copy_particles(p0, sp.N);
    Particle *p_cell = copy_particles(p0, sp.N);
    Particle *p_nbl  = copy_particles(p0, sp.N);
    Particle *p_man  = copy_particles(p0, sp.N);
    if (!p_full || !p_cell || !p_nbl || !p_man) {
        fprintf(stderr, "Allocation failure\n");
        free(p0); free(p_full); free(p_cell); free(p_nbl); free(p_man); return 1;
    }

    // ------------------------------------------------
    // AUTOTUNING
    // ------------------------------------------------
    printf("\n================ BEGINNING AUTOTUNING ================\n");

    printf("\n================ FULL O(N^2) GPU (%d steps) ================\n",
           AUTOTUNE_N_TIMESTEPS);
    FILE *f_full = open_energy_csv("output/full_energies.csv");
    if (!f_full) { free(p0); free(p_full); free(p_cell); free(p_nbl); free(p_man); return 1; }
    double time_full = cuda_full(&sp, p_full, f_full, AUTOTUNE_N_TIMESTEPS);
    fclose(f_full);
    io_write_pdb("output/full_positions.pdb", p_full, sp.N);
    printf("Full O(N^2) done.   Time: %.6f s\n", time_full);

    printf("\n================ CELL-LIST GPU (%d steps) ================\n",
           AUTOTUNE_N_TIMESTEPS);
    FILE *f_cell = open_energy_csv("output/cell_energies.csv");
    if (!f_cell) { free(p0); free(p_full); free(p_cell); free(p_nbl); free(p_man); return 1; }
    double time_cell = cuda_cell_list(&sp, p_cell, f_cell, AUTOTUNE_N_TIMESTEPS);
    fclose(f_cell);
    io_write_pdb("output/cell_positions.pdb", p_cell, sp.N);
    printf("Cell-list done.     Time: %.6f s\n", time_cell);

    printf("\n================ NEIGHBOR-LIST GPU (%d steps) ================\n",
           AUTOTUNE_N_TIMESTEPS);
    FILE *f_nbl = open_energy_csv("output/nbl_energies.csv");
    if (!f_nbl) { free(p0); free(p_full); free(p_cell); free(p_nbl); free(p_man); return 1; }
    double time_nbl = cuda_nbl(&sp, p_nbl, f_nbl, AUTOTUNE_N_TIMESTEPS);
    fclose(f_nbl);
    io_write_pdb("output/nbl_positions.pdb", p_nbl, sp.N);
    printf("Neighbor-list done. Time: %.6f s\n", time_nbl);

    printf("\n================ MANHATTAN CELL-LIST GPU (%d steps) ================\n",
           AUTOTUNE_N_TIMESTEPS);
    FILE *f_man = open_energy_csv("output/manhattan_energies.csv");
    if (!f_man) { free(p0); free(p_full); free(p_cell); free(p_nbl); free(p_man); return 1; }
    double time_man = cuda_manhattan(&sp, p_man, f_man, AUTOTUNE_N_TIMESTEPS);
    fclose(f_man);
    io_write_pdb("output/manhattan_positions.pdb", p_man, sp.N);
    printf("Manhattan done.     Time: %.6f s\n", time_man);

    // ------------------------------------------------
    // TIMING COMPARISON
    // ------------------------------------------------
    printf("\n========================================\n");
    printf("PERFORMANCE COMPARISON\n");
    printf("========================================\n");
    printf("Full O(N^2):   %.6f s\n", time_full);
    printf("Cell-list:     %.6f s  (%.2fx)\n", time_cell, time_full / time_cell);
    printf("Neighbor-list: %.6f s  (%.2fx)\n", time_nbl,  time_full / time_nbl);
    printf("Manhattan:     %.6f s  (%.2fx)\n", time_man,  time_full / time_man);
    printf("========================================\n");

    // Determine fastest
    int fastest_method = 0;
    double fastest_time = time_full;
    const char *fastest_name = "FULL O(N^2) GPU";

    if (time_cell < fastest_time) {
        fastest_time   = time_cell;
        fastest_method = 1;
        fastest_name   = "CELL-LIST GPU";
    }
    if (time_nbl < fastest_time) {
        fastest_time   = time_nbl;
        fastest_method = 2;
        fastest_name   = "NEIGHBOR-LIST GPU";
    }
    if (time_man < fastest_time) {
        fastest_time   = time_man;
        fastest_method = 3;
        fastest_name   = "MANHATTAN GPU";
    }
    printf("FASTEST: %s\n", fastest_name);

    // ------------------------------------------------
    // FULL SIMULATION with fastest algorithm
    // ------------------------------------------------
    printf("\n================ FULL RUN — %s (%d steps) ================\n",
           fastest_name, USER_N_TIMESTEPS);

    Particle *p_final = copy_particles(p0, sp.N);
    if (!p_final) { free(p0); free(p_full); free(p_cell); free(p_nbl); free(p_man); return 1; }

    FILE *f_final = open_energy_csv("output/final_energies.csv");
    if (!f_final) { free(p0); free(p_full); free(p_cell); free(p_nbl); free(p_man); free(p_final); return 1; }

    double time_final;
    if (fastest_method == 0)
        time_final = cuda_full(&sp, p_final, f_final, USER_N_TIMESTEPS);
    else if (fastest_method == 1)
        time_final = cuda_cell_list(&sp, p_final, f_final, USER_N_TIMESTEPS);
    else if (fastest_method == 2)
        time_final = cuda_nbl(&sp, p_final, f_final, USER_N_TIMESTEPS);
    else
        time_final = cuda_manhattan(&sp, p_final, f_final, USER_N_TIMESTEPS);

    fclose(f_final);
    io_write_pdb("output/final_positions.pdb", p_final, sp.N);
    printf("Full run done. Time: %.6f s\n", time_final);
    printf("Output: output/final_energies.csv, output/final_positions.pdb\n");

    free(p0); free(p_full); free(p_cell); free(p_nbl); free(p_man); free(p_final);
    printf("\nAll simulations complete.\n");
    return 0;
}
