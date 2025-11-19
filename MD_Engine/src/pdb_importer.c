#include "pdb_importer.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

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
