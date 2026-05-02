#pragma once

#include <cuda_runtime.h>
#include <cstdio>

// ---------------------------------------------------------------------------
// Error checking
// ---------------------------------------------------------------------------

#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = (call); \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                    __FILE__, __LINE__, cudaGetErrorString(err)); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

// ---------------------------------------------------------------------------
// Cell boundary finder
// Requires cellIDs to be sorted.  Writes the half-open interval
// [cellStart[c], cellEnd[c]) for each occupied cell c.
// Caller should cudaMemset cellStart/cellEnd to -1 before launching.
// static so each translation unit gets its own copy (safe with -rdc=true).
// ---------------------------------------------------------------------------

static __global__ void findCellBounds(const int *cellIDs,
                                      int *cellStart, int *cellEnd, int N) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    int current = cellIDs[i];
    if (i == 0) {
        cellStart[current] = 0;
    } else {
        int prev = cellIDs[i - 1];
        if (current != prev) {
            cellEnd[prev]      = i;
            cellStart[current] = i;
        }
    }
    if (i == N - 1) cellEnd[current] = N;
}
