#include "md_io.h"
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// ---------------------------------------------------------------------------
// INPUT
// ---------------------------------------------------------------------------

static int binary_importer(const char *path,
                            Particle *p, int n_atoms,
                            double *xmin, double *xmax,
                            double *ymin, double *ymax,
                            double *zmin, double *zmax)
{
    FILE *f = fopen(path, "rb");
    if (!f) return -1;

    double *buffer = (double *)malloc(n_atoms * 3 * sizeof(double));
    if (!buffer) { fclose(f); return -1; }

    fread(buffer, sizeof(double), n_atoms * 3, f);

    *xmin = *ymin = *zmin =  1e99;
    *xmax = *ymax = *zmax = -1e99;

    for (int i = 0; i < n_atoms; i++) {
        double x = buffer[i * 3 + 0];
        double y = buffer[i * 3 + 1];
        double z = buffer[i * 3 + 2];

        p[i].x = x; p[i].y = y; p[i].z = z;
        p[i].vx = p[i].vy = p[i].vz = 0.0;
        p[i].fx = p[i].fy = p[i].fz = 0.0;

        if (x < *xmin) *xmin = x; if (x > *xmax) *xmax = x;
        if (y < *ymin) *ymin = y; if (y > *ymax) *ymax = y;
        if (z < *zmin) *zmin = z; if (z > *zmax) *zmax = z;
    }

    free(buffer);
    fclose(f);
    return n_atoms;
}

int read_input(const char *input_path, Particle **p0, SimParams *sp)
{
    printf("Converting input file via Python script...\n");
    char command[256];
    snprintf(command, sizeof(command),
             "./venv/bin/python3 src/py/converter.py %s temp_coords.bin temp_metadata.txt",
             input_path);
    int status = system(command);
    if (status != 0) {
        fprintf(stderr, "Error: Python conversion failed.\n");
        return -1;
    }

    FILE *f_meta = fopen("temp_metadata.txt", "r");
    if (!f_meta) {
        fprintf(stderr, "Metadata not found!\n");
        return -1;
    }
    int N;
    fscanf(f_meta, "%d", &N);
    fclose(f_meta);

    *p0 = (Particle *)malloc(N * sizeof(Particle));

    double xmin, xmax, ymin, ymax, zmin, zmax;
    binary_importer("temp_coords.bin", *p0, N,
                    &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

    system("rm -f temp_coords.bin temp_metadata.txt");

    printf("Loaded %d atoms from %s\n", N, input_path);
    printf("Bounds: dx = %.3f  dy = %.3f  dz = %.3f\n",
           xmax - xmin, ymax - ymin, zmax - zmin);

    sp->N = (size_t)N;
    sp->L = fmax(fmax(xmax - xmin, ymax - ymin), zmax - zmin);
    if (sp->L <= 0.0) sp->L = CONF_BOX_MAX;

    return 0;
}

// ---------------------------------------------------------------------------
// OUTPUT
// ---------------------------------------------------------------------------

void io_write_pdb(const char *path, Particle *p, size_t N)
{
    FILE *f = fopen(path, "w");
    if (!f) { perror("fopen"); return; }

    for (size_t i = 0; i < N; i++) {
        fprintf(f,
            "ATOM  %5zu  X   XXX A   1    %8.3f %8.3f %8.3f  1.00  0.00           X\n",
            i + 1, p[i].x, p[i].y, p[i].z);
    }

    fprintf(f, "END\n");
    fclose(f);
}
