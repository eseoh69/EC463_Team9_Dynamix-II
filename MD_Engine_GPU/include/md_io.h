#ifndef MD_IO_H
#define MD_IO_H

#ifdef __cplusplus
extern "C" {
#endif

#include "md.h"

// parse PDB via MDAnalysis converter, allocates *p0 and fills *sp (N and L)
int read_input(const char *input_path, Particle **p0, SimParams *sp);

// write particle positions to PDB file
void io_write_pdb(const char *path, Particle *p, size_t N);

#ifdef __cplusplus
}
#endif

#endif
