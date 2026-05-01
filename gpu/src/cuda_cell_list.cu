#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <thrust/sequence.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <chrono>
#include <cstdio>

#include "md.h"
#include "md_io.h"
#include "cuda_cell_list.h"
#include "config.h"
#include "cuda_common.cuh"

// ---------------------------------------------------------------------------
// Device helpers
// ---------------------------------------------------------------------------

// Helper for 3D vector subtraction + Minimum Image Convention
__device__ double3 getDiffPBC(double3 p1, double3 p2, double simSize) {
    double3 diff = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
    if (diff.x >  simSize * 0.5) diff.x -= simSize;
    if (diff.x < -simSize * 0.5) diff.x += simSize;
    if (diff.y >  simSize * 0.5) diff.y -= simSize;
    if (diff.y < -simSize * 0.5) diff.y += simSize;
    if (diff.z >  simSize * 0.5) diff.z -= simSize;
    if (diff.z < -simSize * 0.5) diff.z += simSize;
    return diff;
}

__device__ double3 computeLJForce(double3 p1, double3 p2, double simSize,
                                   double sigma, double epsilon, double rc2) {
    double3 diff = getDiffPBC(p1, p2, simSize);
    double distSq = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;

    // Cutoff check
    if (distSq < 1e-30 || distSq > rc2) return {0, 0, 0};

    double inv_r2   = 1.0 / distSq;
    double sig2_r2  = (sigma * sigma) * inv_r2;
    double sr6      = sig2_r2 * sig2_r2 * sig2_r2;
    double sr12     = sr6 * sr6;
    double forceMag = 24.0 * epsilon * (2.0*sr12 - sr6) * inv_r2;

    return {-diff.x * forceMag, -diff.y * forceMag, -diff.z * forceMag};
}

// ---------------------------------------------------------------------------
// Kernels
// ---------------------------------------------------------------------------

// Integration
__global__ void verletStage1(Particle_cuda *p, double dt, double simSize, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    double invM = 1.0; // assume mass 1
    p[i].pos.x += p[i].vel.x * dt + 0.5 * p[i].force.x * invM * dt * dt;
    p[i].pos.y += p[i].vel.y * dt + 0.5 * p[i].force.y * invM * dt * dt;
    p[i].pos.z += p[i].vel.z * dt + 0.5 * p[i].force.z * invM * dt * dt;

    p[i].vel.x += 0.5 * p[i].force.x * invM * dt;
    p[i].vel.y += 0.5 * p[i].force.y * invM * dt;
    p[i].vel.z += 0.5 * p[i].force.z * invM * dt;

    p[i].pos.x = fmod(p[i].pos.x, simSize); if (p[i].pos.x < 0) p[i].pos.x += simSize;
    p[i].pos.y = fmod(p[i].pos.y, simSize); if (p[i].pos.y < 0) p[i].pos.y += simSize;
    p[i].pos.z = fmod(p[i].pos.z, simSize); if (p[i].pos.z < 0) p[i].pos.z += simSize;
}

__global__ void cellComputeKE(const Particle_cuda *p, double *ke_out, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    ke_out[i] = 0.5 * (p[i].vel.x*p[i].vel.x + p[i].vel.y*p[i].vel.y + p[i].vel.z*p[i].vel.z);
}

__global__ void verletStage2(Particle_cuda *p, double dt, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    double invM = 1.0;
    p[i].vel.x += 0.5 * p[i].force.x * invM * dt;
    p[i].vel.y += 0.5 * p[i].force.y * invM * dt;
    p[i].vel.z += 0.5 * p[i].force.z * invM * dt;
}

// 2. CELL ID CALCULATION
__global__ void calculateCellIDs(Particle_cuda *p, int *cellIDs, int3 gridDim,
                                   double cellSize, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    int cx = (int)(p[i].pos.x / cellSize);
    int cy = (int)(p[i].pos.y / cellSize);
    int cz = (int)(p[i].pos.z / cellSize);

    // Wrap cell indices to handle positions outside [0, L]
    cx = ((cx % gridDim.x) + gridDim.x) % gridDim.x;
    cy = ((cy % gridDim.y) + gridDim.y) % gridDim.y;
    cz = ((cz % gridDim.z) + gridDim.z) % gridDim.z;

    cellIDs[i] = (cz * gridDim.y * gridDim.x) + (cy * gridDim.x) + cx;
}

// 3. FORCE KERNEL (3D Neighbor Search)
__global__ void computeForces3D(Particle_cuda *p, int *cellStart, int *cellEnd,
                                  int *d_indices, int3 gridDim,
                                  double simSize, double cellSize, int n,
                                  double sigma, double epsilon, double rc2) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    int    real_idx = d_indices[i];
    double3 myPos   = p[real_idx].pos;
    double3 totalF  = {0, 0, 0};

    int cx = (int)(myPos.x / cellSize);
    int cy = (int)(myPos.y / cellSize);
    int cz = (int)(myPos.z / cellSize);

    for (int sz = cz - 1; sz <= cz + 1; ++sz) {
        for (int sy = cy - 1; sy <= cy + 1; ++sy) {
            for (int sx = cx - 1; sx <= cx + 1; ++sx) {
                // Handling PBC for neighbor cells
                int nz = (sz + gridDim.z) % gridDim.z;
                int ny = (sy + gridDim.y) % gridDim.y;
                int nx = (sx + gridDim.x) % gridDim.x;

                int neighborCell = (nz * gridDim.y * gridDim.x) + (ny * gridDim.x) + nx;
                int start = cellStart[neighborCell];
                int end   = cellEnd[neighborCell];

                if (start == -1) continue;

                for (int j_ptr = start; j_ptr < end; ++j_ptr) {
                    if (i == j_ptr) continue;
                    // Look up the neighbor's REAL memory location
                    int neighbor_real_idx = d_indices[j_ptr];
                    double3 f = computeLJForce(myPos, p[neighbor_real_idx].pos,
                                               simSize, sigma, epsilon, rc2);
                    totalF.x += f.x;
                    totalF.y += f.y;
                    totalF.z += f.z;
                }
            }
        }
    }

    p[real_idx].force = totalF; // Write back to the original memory spot
}

// Potential energy per particle (each pair counted once via i < j in sorted order)
__global__ void computePE(Particle_cuda *p, int *cellStart, int *cellEnd,
                           int *d_indices, int3 gridDim,
                           double simSize, double cellSize, int n,
                           double sigma, double epsilon, double rc2,
                           double *pe_out) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    int     real_idx = d_indices[i];
    double3 myPos    = p[real_idx].pos;
    double  pe       = 0.0;

    int cx = (int)(myPos.x / cellSize);
    int cy = (int)(myPos.y / cellSize);
    int cz = (int)(myPos.z / cellSize);

    for (int sz = cz - 1; sz <= cz + 1; ++sz) {
        for (int sy = cy - 1; sy <= cy + 1; ++sy) {
            for (int sx = cx - 1; sx <= cx + 1; ++sx) {
                int nz = (sz + gridDim.z) % gridDim.z;
                int ny = (sy + gridDim.y) % gridDim.y;
                int nx = (sx + gridDim.x) % gridDim.x;

                int neighborCell = (nz * gridDim.y * gridDim.x) + (ny * gridDim.x) + nx;
                int start = cellStart[neighborCell];
                int end   = cellEnd[neighborCell];

                if (start == -1) continue;

                for (int j_ptr = start; j_ptr < end; ++j_ptr) {
                    if (j_ptr <= i) continue;  // count each pair once
                    int neighbor_real_idx = d_indices[j_ptr];
                    double3 diff = getDiffPBC(myPos, p[neighbor_real_idx].pos, simSize);
                    double  r2   = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
                    if (r2 < 1e-30 || r2 > rc2) continue;
                    double inv_r2  = 1.0 / r2;
                    double sig2_r2 = sigma * sigma * inv_r2;
                    double sr6     = sig2_r2 * sig2_r2 * sig2_r2;
                    double sr12    = sr6 * sr6;
                    pe += 4.0 * epsilon * (sr12 - sr6);
                }
            }
        }
    }

    pe_out[i] = pe;
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

double cuda_cell_list(const SimParams *sp, Particle *p, FILE *f_energy, int nsteps) {
    const int    numParticles = (int)sp->N;
    const double simSize      = sp->L;
    const double cellSize     = sp->rc;
    int cellDim = (int)(simSize / cellSize);
    if (cellDim < 1) cellDim = 1;  // guard against L < rc

    const int3   gridDim      = {cellDim, cellDim, cellDim};
    const int    numCells     = cellDim * cellDim * cellDim;

    // Particle and Particle_cuda share the same memory layout (9 doubles)
    Particle_cuda *h_particles = reinterpret_cast<Particle_cuda *>(p);

    // Allocate device memory
    Particle_cuda *d_particles;
    int    *d_cellIDs, *d_cellStart, *d_cellEnd, *d_indices;
    double *d_pe_out;

    CUDA_CHECK(cudaMalloc(&d_particles, numParticles * sizeof(Particle_cuda)));
    CUDA_CHECK(cudaMalloc(&d_cellIDs,   numParticles * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_cellStart, numCells     * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_cellEnd,   numCells     * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_indices,   numParticles * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_pe_out,    numParticles * sizeof(double)));

    thrust::sequence(thrust::device, d_indices, d_indices + numParticles);
    CUDA_CHECK(cudaMemcpy(d_particles, h_particles, numParticles * sizeof(Particle_cuda),
               cudaMemcpyHostToDevice));

    const int blockSize = GPU_BLOCK_SIZE;
    const int numBlocks = (numParticles + blockSize - 1) / blockSize;

    // Compute initial forces before first integration step
    {
        thrust::device_ptr<int> d_keys(d_cellIDs);
        thrust::device_ptr<int> d_vals(d_indices);

        calculateCellIDs<<<numBlocks, blockSize>>>(d_particles, d_cellIDs, gridDim,
                                                    cellSize, numParticles);
        CUDA_CHECK(cudaDeviceSynchronize());

        thrust::sequence(thrust::device, d_indices, d_indices + numParticles);
        thrust::sort_by_key(d_keys, d_keys + numParticles, d_vals);

        CUDA_CHECK(cudaMemset(d_cellStart, -1, numCells * sizeof(int)));
        CUDA_CHECK(cudaMemset(d_cellEnd,   -1, numCells * sizeof(int)));

        findCellBounds<<<numBlocks, blockSize>>>(d_cellIDs, d_cellStart, d_cellEnd, numParticles);
        CUDA_CHECK(cudaDeviceSynchronize());

        computeForces3D<<<numBlocks, blockSize>>>(d_particles, d_cellStart, d_cellEnd,
                                                   d_indices, gridDim, simSize, cellSize,
                                                   numParticles, sp->sigma, sp->epsilon, sp->rc2);
        CUDA_CHECK(cudaDeviceSynchronize());
    }

    auto start = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < nsteps; step++) {
        verletStage1<<<numBlocks, blockSize>>>(d_particles, sp->dt, simSize, numParticles);
        CUDA_CHECK(cudaDeviceSynchronize());

        thrust::device_ptr<int> d_keys(d_cellIDs);
        thrust::device_ptr<int> d_vals(d_indices);

        calculateCellIDs<<<numBlocks, blockSize>>>(d_particles, d_cellIDs, gridDim,
                                                    cellSize, numParticles);
        CUDA_CHECK(cudaDeviceSynchronize());

        thrust::sequence(thrust::device, d_indices, d_indices + numParticles);
        thrust::sort_by_key(d_keys, d_keys + numParticles, d_vals);

        CUDA_CHECK(cudaMemset(d_cellStart, -1, numCells * sizeof(int)));
        CUDA_CHECK(cudaMemset(d_cellEnd,   -1, numCells * sizeof(int)));
        findCellBounds<<<numBlocks, blockSize>>>(d_cellIDs, d_cellStart, d_cellEnd, numParticles);
        CUDA_CHECK(cudaDeviceSynchronize());

        computeForces3D<<<numBlocks, blockSize>>>(d_particles, d_cellStart, d_cellEnd,
                                                   d_indices, gridDim, simSize, cellSize,
                                                   numParticles, sp->sigma, sp->epsilon, sp->rc2);
        CUDA_CHECK(cudaDeviceSynchronize());

        verletStage2<<<numBlocks, blockSize>>>(d_particles, sp->dt, numParticles);
        CUDA_CHECK(cudaDeviceSynchronize());

        if (f_energy) {
            cellComputeKE<<<numBlocks, blockSize>>>(d_particles, d_pe_out, numParticles);
            double K = thrust::reduce(thrust::device_ptr<double>(d_pe_out),
                                      thrust::device_ptr<double>(d_pe_out + numParticles),
                                      0.0, thrust::plus<double>());

            computePE<<<numBlocks, blockSize>>>(d_particles, d_cellStart, d_cellEnd,
                                                d_indices, gridDim, simSize, cellSize,
                                                numParticles, sp->sigma, sp->epsilon, sp->rc2,
                                                d_pe_out);
            double U = thrust::reduce(thrust::device_ptr<double>(d_pe_out),
                                      thrust::device_ptr<double>(d_pe_out + numParticles),
                                      0.0, thrust::plus<double>());

            fprintf(f_energy, "%d,%.10f,%.10f,%.10f\n", step, K, U, K + U);
        }
    }

    cudaDeviceSynchronize();
    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = stop - start;

    // Copy results back
    cudaMemcpy(h_particles, d_particles, numParticles * sizeof(Particle_cuda),
               cudaMemcpyDeviceToHost);

    // Cleanup
    cudaFree(d_particles);
    cudaFree(d_cellIDs);
    cudaFree(d_cellStart);
    cudaFree(d_cellEnd);
    cudaFree(d_indices);
    cudaFree(d_pe_out);

    return elapsed.count();
}
