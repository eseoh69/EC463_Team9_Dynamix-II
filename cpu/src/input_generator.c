/*
 * generate_particles.c
 *
 * Generates random-particle PDB input files for MD simulations.
 * Scales the cubic box length so that number density remains constant:
 *
 *   rho = N / L^3  = const  =>  L(N) = L0 * (N / N0)^(1/3)
 *
 * Reference point (from EC463_Team9_Dynamix-II input):
 *   N0 = 1024  particles,  L0 = 300.0 Angstroms
 *
 * Particle counts generated: 1024, 2048, 4096, 8192, 16384, 32768
 *
 * Usage:
 *   gcc -O2 -o generate_particles generate_particles.c -lm
 *   ./generate_particles
 *
 * Output files are written to the current directory:
 *   random_particles-<N>.pdb
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
 
/* ── tuneable constants ────────────────────────────────────────────────── */
#define N0      1024        /* reference particle count                     */
#define L0      300.0       /* reference box length (Angstroms)             */
#define OCCUPANCY 1.00
#define BFACTOR   0.00
/* ───────────────────────────────────────────────────────────────────────── */
 
/* Uniform random double in [0, range) */
static inline double rand_coord(double range)
{
    return range * ((double)rand() / ((double)RAND_MAX + 1.0));
}
 
/*
 * generate_pdb  –  writes one PDB file for N particles in a cube of
 *                  side length box_len.
 */
static int generate_pdb(int N, double box_len)
{
    char filename[64];
    snprintf(filename, sizeof(filename), "random_particles-%d.pdb", N);
 
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open %s for writing\n", filename);
        return -1;
    }
 
    /*
     * CRYST1 record – tells MD engines the periodic box dimensions.
     * Angles are 90.0 for a rectangular (cubic) box.
     */
    fprintf(fp,
            "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
            box_len, box_len, box_len, 90.0, 90.0, 90.0);
 
    for (int i = 1; i <= N; i++) {
        double x = rand_coord(box_len);
        double y = rand_coord(box_len);
        double z = rand_coord(box_len);
 
        /*
         * Standard PDB ATOM record (columns strictly fixed-width):
         *
         * cols  1- 6  record type  "ATOM  "
         * cols  7-11  serial       right-justified integer
         * col  12     space
         * cols 13-16  atom name    " X  "
         * col  17     alt loc      " "
         * cols 18-20  res name     "XXX"
         * col  21     space
         * col  22     chain ID     "A"
         * cols 23-26  res seq      right-justified integer
         * col  27     i-code       " "
         * cols 28-30  spaces
         * cols 31-38  x            %8.3f
         * cols 39-46  y            %8.3f
         * cols 47-54  z            %8.3f
         * cols 55-60  occupancy    %6.2f
         * cols 61-66  temp factor  %6.2f
         * cols 77-78  element      " X"
         */
        fprintf(fp,
                "ATOM  %5d  X   XXX A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f"
                "          %2s\n",
                i, 1,
                x, y, z,
                OCCUPANCY, BFACTOR,
                "X");
    }
 
    fprintf(fp, "END\n");
    fclose(fp);
    return 0;
}
 
int main(void)
{
    /* Particle counts to generate (powers of 2 from 1024 to 32768) */
    int counts[] = { 1024, 2048, 4096, 8192, 16384, 32768 };
    int n_files  = (int)(sizeof(counts) / sizeof(counts[0]));
 
    srand((unsigned int)time(NULL));   /* seed RNG from wall clock */
 
    printf("Generating PDB files (reference: N0=%d, L0=%.1f Å)\n\n",
           N0, L0);
    printf("  %-10s  %-14s  %-20s\n",
           "N", "Box length (Å)", "Output file");
    printf("  %-10s  %-14s  %-20s\n",
           "----------", "--------------", "--------------------");
 
    for (int f = 0; f < n_files; f++) {
        int    N       = counts[f];
        double box_len = L0 * cbrt((double)N / (double)N0);
 
        printf("  %-10d  %-14.3f  random_particles-%d.pdb\n",
               N, box_len, N);
 
        if (generate_pdb(N, box_len) != 0) {
            fprintf(stderr, "Failed to generate file for N=%d\n", N);
            return EXIT_FAILURE;
        }
    }
 
    printf("\nDone. %d files written.\n", n_files);
    return EXIT_SUCCESS;
}
 