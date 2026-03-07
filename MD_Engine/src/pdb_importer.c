#include "pdb_importer.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "config.h"
#include <math.h>

int pdb_importer(const char *path,
                 Particle *p, int maxN,
                 double *xmin, double *xmax,
                 double *ymin, double *ymax,
                 double *zmin, double *zmax)
{
    FILE *f = fopen(path, "r");
    if (!f) return -1;

    char line[256];
    int count = 0;

    *xmin = *ymin = *zmin = 1e99;
    *xmax = *ymax = *zmax = -1e99;

    while (fgets(line, sizeof(line), f) && count < maxN) {
        if (!strncmp(line, "ATOM", 4)) {

            double x = atof(line + 30);
            double y = atof(line + 38);
            double z = atof(line + 46);

            p[count].x = x;
            p[count].y = y;
            p[count].z = z;

            p[count].vx = p[count].vy = p[count].vz = 0;
            p[count].fx = p[count].fy = p[count].fz = 0;

            if (x < *xmin) *xmin = x;
            if (x > *xmax) *xmax = x;
            if (y < *ymin) *ymin = y;
            if (y > *ymax) *ymax = y;
            if (z < *zmin) *zmin = z;
            if (z > *zmax) *zmax = z;

            count++;
        }
    }

    fclose(f);
    return count;
}

int binary_importer(const char *path, 
                    Particle *p, int n_atoms,
                    double *xmin, double *xmax, 
                    double *ymin, double *ymax, 
                    double *zmin, double *zmax) 
{
    FILE *f = fopen(path, "rb");
    if (!f) return -1;

    // We need a temporary buffer to read the (x,y,z) triplets 
    // because the 'Particle' struct has extra fields (v and f) 
    // that aren't in the binary file.
    double *buffer = (double *)malloc(n_atoms * 3 * sizeof(double));
    if (!buffer) { fclose(f); return -1; }

    // Read all coordinates in one shot
    fread(buffer, sizeof(double), n_atoms * 3, f);

    *xmin = *ymin = *zmin = 1e99;
    *xmax = *ymax = *zmax = -1e99;

    for (int i = 0; i < n_atoms; i++) {
        double x = buffer[i * 3 + 0];
        double y = buffer[i * 3 + 1];
        double z = buffer[i * 3 + 2];

        p[i].x = x;
        p[i].y = y;
        p[i].z = z;

        // Initialize physics vectors
        p[i].vx = p[i].vy = p[i].vz = 0;
        p[i].fx = p[i].fy = p[i].fz = 0;

        // Update bounds
        if (x < *xmin) *xmin = x; if (x > *xmax) *xmax = x;
        if (y < *ymin) *ymin = y; if (y > *ymax) *ymax = y;
        if (z < *zmin) *zmin = z; if (z > *zmax) *zmax = z;
    }

    free(buffer);
    fclose(f);
    return n_atoms;
}

int read_input (const char * input_path,Particle ** p0, SimParams * sp){

   printf("Converting input file via Python script...\n");
    char command[256];
    sprintf(command, "./venv/bin/python3 src/py/converter.py %s temp_coords.bin temp_metadata.txt", input_path);
    int status = system(command);

    if (status != 0) {
        fprintf(stderr, "Error: Python conversion failed.\n");
        return -1;
    }

    FILE *f_meta = fopen("temp_metadata.txt", "r");
    int N;
    if (f_meta) {
        fscanf(f_meta, "%d", &N);
        fclose(f_meta);
    } else {
        fprintf(stderr, "Metadata not found!\n");
        return -1;
    }

    // Now you can allocate exactly what you need
    *p0 = malloc(N * sizeof(Particle));

    // Now you know temp_coords.bin exists and is ready
    int foo;

    double xmin, xmax, ymin, ymax, zmin, zmax;
    foo = binary_importer("temp_coords.bin", *p0, N,
                         &xmin, &xmax,
                         &ymin, &ymax,
                         &zmin, &zmax);
    (void)foo; // Suppress unused variable warning
    
    if (N <= 0) {
        fprintf(stderr, "Failed to read PDB: %s\n", input_path);
        free(*p0);
        return -1;
    }

    printf("deleting temp files \n");
    status = system("rm -rf temp_coords.bin");
    status = system("rm -rf temp_metadata.txt");

    printf("Loaded %d atoms from %s\n", N, input_path);
    printf("Bounds: dx = %.3f  dy = %.3f  dz = %.3f\n",
           xmax - xmin, ymax - ymin, zmax - zmin);

    // set N and L
    sp->L = fmax(fmax(xmax - xmin, ymax - ymin), zmax - zmin);
    if (sp->L <= 0.0) sp->L = CONF_BOX_MAX;

    sp->N = (size_t)N;
    return 0;
    
}