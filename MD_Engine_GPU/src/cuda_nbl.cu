#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/sequence.h>
#include <thrust/reduce.h>
#include <chrono>
#include <cstdio>

#include "md.h"
#include "md_io.h"
#include "cuda_nbl.h"
#include "config.h"
#include "cuda_common.cuh"

// ---------------------------------------------------------------------------
// Kernels
// ---------------------------------------------------------------------------

__global__ void calcCellIDs(
    const double *x, const double *y, const double *z,
    int *cellIDs, int nc, double cellSize, int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    int cx = (int)(x[i] / cellSize);
    int cy = (int)(y[i] / cellSize);
    int cz = (int)(z[i] / cellSize);

    // clamp to valid range
    if (cx >= nc) cx = nc - 1;  if (cx < 0) cx = 0;
    if (cy >= nc) cy = nc - 1;  if (cy < 0) cy = 0;
    if (cz >= nc) cz = nc - 1;  if (cz < 0) cz = 0;

    cellIDs[i] = cx + nc * (cy + nc * cz);
}

__global__ void buildNeighborList(
    const double *x, const double *y, const double *z,
    const int *cellStart, const int *cellEnd,
    const int *sortedIndices,
    int *nb, int *nbCount,
    int nc, double cellSize, double L, double rc_skin2,
    int N, int maxNeighbors)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    double xi = x[i], yi = y[i], zi = z[i];
    double halfL = 0.5 * L;

    int cx = (int)(xi / cellSize); if (cx >= nc) cx = nc-1; if (cx < 0) cx = 0;
    int cy = (int)(yi / cellSize); if (cy >= nc) cy = nc-1; if (cy < 0) cy = 0;
    int cz = (int)(zi / cellSize); if (cz >= nc) cz = nc-1; if (cz < 0) cz = 0;

    int count = 0;
    int base  = i * maxNeighbors;

    for (int dz = -1; dz <= 1; dz++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                int nx = (cx + dx + nc) % nc;
                int ny = (cy + dy + nc) % nc;
                int nz = (cz + dz + nc) % nc;

                int cellIdx = nx + nc * (ny + nc * nz);
                int start   = cellStart[cellIdx];
                int end     = cellEnd[cellIdx];
                if (start == -1) continue;

                for (int k = start; k < end; k++) {
                    int j = sortedIndices[k];
                    if (j <= i) continue;  // half-shell

                    double dxij = xi - x[j];
                    double dyij = yi - y[j];
                    double dzij = zi - z[j];

                    if (dxij >  halfL) dxij -= L; else if (dxij < -halfL) dxij += L;
                    if (dyij >  halfL) dyij -= L; else if (dyij < -halfL) dyij += L;
                    if (dzij >  halfL) dzij -= L; else if (dzij < -halfL) dzij += L;

                    double r2 = dxij*dxij + dyij*dyij + dzij*dzij;
                    if (r2 <= rc_skin2 && count < maxNeighbors) {
                        nb[base + count] = j;
                        count++;
                    }
                }
            }
        }
    }

    nbCount[i] = count;
}

__global__ void ljForceKernel(
    const double *x, const double *y, const double *z,
    double *fx, double *fy, double *fz,
    double *U_per_particle,
    const int *nb, const int *nbCount,
    int N, int maxNeighbors,
    double eps, double sig, double L, double rc2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    double sig2  = sig * sig;
    double halfL = 0.5 * L;

    double xi = x[i], yi = y[i], zi = z[i];
    double fx_i = 0.0, fy_i = 0.0, fz_i = 0.0, U_i = 0.0;

    int base         = i * maxNeighbors;
    int numNeighbors = nbCount[i];

    for (int k = 0; k < numNeighbors; k++) {
        int j = nb[base + k];

        double dx = xi - x[j];
        double dy = yi - y[j];
        double dz = zi - z[j];

        if (dx >  halfL) dx -= L; else if (dx < -halfL) dx += L;
        if (dy >  halfL) dy -= L; else if (dy < -halfL) dy += L;
        if (dz >  halfL) dz -= L; else if (dz < -halfL) dz += L;

        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1e-30 || r2 > rc2) continue;

        double inv_r2  = 1.0 / r2;
        double sig2_r2 = sig2 * inv_r2;
        double sr6     = sig2_r2 * sig2_r2 * sig2_r2;
        double sr12    = sr6 * sr6;
        double F_over_r = 24.0 * eps * (2.0*sr12 - sr6) * inv_r2;

        double fxij = F_over_r * dx;
        double fyij = F_over_r * dy;
        double fzij = F_over_r * dz;

        // accumulate force on i locally
        fx_i += fxij;
        fy_i += fyij;
        fz_i += fzij;

        // Newton's 3rd: force on j via atomicAdd (SM 6.0+ supports native double atomicAdd)
        atomicAdd(&fx[j], -fxij);
        atomicAdd(&fy[j], -fyij);
        atomicAdd(&fz[j], -fzij);

        U_i += 4.0 * eps * (sr12 - sr6);
    }

    atomicAdd(&fx[i], fx_i);
    atomicAdd(&fy[i], fy_i);
    atomicAdd(&fz[i], fz_i);

    U_per_particle[i] = U_i;
}

__global__ void nbl_verletStage1(
    double *x, double *y, double *z,
    double *vx, double *vy, double *vz,
    const double *fx, const double *fy, const double *fz,
    double dt, double L, int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    double half_dt = 0.5 * dt;

    vx[i] += half_dt * fx[i];
    vy[i] += half_dt * fy[i];
    vz[i] += half_dt * fz[i];

    x[i] += dt * vx[i];
    y[i] += dt * vy[i];
    z[i] += dt * vz[i];

    // fmod-based PBC wrap — handles multi-box jumps from large forces
    x[i] = fmod(x[i], L); if (x[i] < 0.0) x[i] += L;
    y[i] = fmod(y[i], L); if (y[i] < 0.0) y[i] += L;
    z[i] = fmod(z[i], L); if (z[i] < 0.0) z[i] += L;
}

__global__ void nbl_verletStage2(
    double *vx, double *vy, double *vz,
    const double *fx, const double *fy, const double *fz,
    double dt, int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    double half_dt = 0.5 * dt;
    vx[i] += half_dt * fx[i];
    vy[i] += half_dt * fy[i];
    vz[i] += half_dt * fz[i];
}

__global__ void computeKE(
    const double *vx, const double *vy, const double *vz,
    double *K_per_particle, int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    K_per_particle[i] = 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
}

__global__ void computeDisplacement(
    const double *x, const double *y, const double *z,
    const double *px, const double *py, const double *pz,
    double *disp2, double L, int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    double dx = x[i] - px[i];
    double dy = y[i] - py[i];
    double dz = z[i] - pz[i];
    double halfL = 0.5 * L;

    if (dx >  halfL) dx -= L; else if (dx < -halfL) dx += L;
    if (dy >  halfL) dy -= L; else if (dy < -halfL) dy += L;
    if (dz >  halfL) dz -= L; else if (dz < -halfL) dz += L;

    disp2[i] = dx*dx + dy*dy + dz*dz;
}

__global__ void savePositions(
    const double *x, const double *y, const double *z,
    double *px, double *py, double *pz, int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    px[i] = x[i]; py[i] = y[i]; pz[i] = z[i];
}

// ---------------------------------------------------------------------------
// Host helpers
// ---------------------------------------------------------------------------

static void buildCellStructure(
    const double *d_x, const double *d_y, const double *d_z,
    int *d_cellIDs, int *d_sortedIdx, int *d_cellStart, int *d_cellEnd,
    int N, int nc, int numCells, double cellSize, int numBlocks)
{
    calcCellIDs<<<numBlocks, GPU_BLOCK_SIZE>>>(
        d_x, d_y, d_z, d_cellIDs, nc, cellSize, N);

    thrust::sequence(thrust::device, d_sortedIdx, d_sortedIdx + N);
    thrust::sort_by_key(thrust::device_ptr<int>(d_cellIDs),
                        thrust::device_ptr<int>(d_cellIDs + N),
                        thrust::device_ptr<int>(d_sortedIdx));

    CUDA_CHECK(cudaMemset(d_cellStart, -1, numCells * sizeof(int)));
    CUDA_CHECK(cudaMemset(d_cellEnd,   -1, numCells * sizeof(int)));
    findCellBounds<<<numBlocks, GPU_BLOCK_SIZE>>>(
        d_cellIDs, d_cellStart, d_cellEnd, N);
}

static int needsRebuild(
    const double *d_x, const double *d_y, const double *d_z,
    const double *d_px, const double *d_py, const double *d_pz,
    double *d_disp2, double L, double skin, int N, int numBlocks)
{
    computeDisplacement<<<numBlocks, GPU_BLOCK_SIZE>>>(
        d_x, d_y, d_z, d_px, d_py, d_pz, d_disp2, L, N);

    double maxDisp2 = thrust::reduce(
        thrust::device_ptr<double>(d_disp2),
        thrust::device_ptr<double>(d_disp2 + N),
        0.0, thrust::maximum<double>());

    double threshold = 0.5 * skin;
    return (maxDisp2 > threshold * threshold) ? 1 : 0;
}

static void rebuildNeighborList(
    const double *d_x, const double *d_y, const double *d_z,
    int *d_cellIDs, int *d_sortedIdx, int *d_cellStart, int *d_cellEnd,
    int *d_nb, int *d_nbCount,
    double *d_px, double *d_py, double *d_pz,
    int N, int nc, int numCells, double cellSize, double L,
    double rc_skin2, int numBlocks)
{
    buildCellStructure(d_x, d_y, d_z,
        d_cellIDs, d_sortedIdx, d_cellStart, d_cellEnd,
        N, nc, numCells, cellSize, numBlocks);

    buildNeighborList<<<numBlocks, GPU_BLOCK_SIZE>>>(
        d_x, d_y, d_z,
        d_cellStart, d_cellEnd, d_sortedIdx,
        d_nb, d_nbCount,
        nc, cellSize, L, rc_skin2, N, MAX_NEIGHBORS);

    savePositions<<<numBlocks, GPU_BLOCK_SIZE>>>(
        d_x, d_y, d_z, d_px, d_py, d_pz, N);
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

double cuda_nbl(const SimParams *sp, Particle *p, FILE *f_energy, int nsteps)
{
    const int    N    = (int)sp->N;
    const double L    = sp->L;
    const double rc   = sp->rc;
    const double skin = NBL_SKIN_FACTOR * rc;
    const double rc_skin  = rc + skin;
    const double rc_skin2 = rc_skin * rc_skin;

    int nc       = (int)(L / rc);
    if (nc < 1) nc = 1;
    double cellSize = L / nc;
    int numCells    = nc * nc * nc;

    const int numBlocks = (N + GPU_BLOCK_SIZE - 1) / GPU_BLOCK_SIZE;

    // ---- Allocate device arrays (SoA layout) ----
    double *d_x, *d_y, *d_z;
    double *d_vx, *d_vy, *d_vz;
    double *d_fx, *d_fy, *d_fz;
    double *d_U, *d_K;
    double *d_px, *d_py, *d_pz, *d_disp2;
    int    *d_cellIDs, *d_sortedIdx, *d_cellStart, *d_cellEnd;
    int    *d_nb, *d_nbCount;

    CUDA_CHECK(cudaMalloc(&d_x,  N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_y,  N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_z,  N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_vx, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_vy, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_vz, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_fx, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_fy, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_fz, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_U,  N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_K,  N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_px, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_py, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_pz, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_disp2, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_cellIDs,  N * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_sortedIdx, N * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_cellStart, numCells * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_cellEnd,   numCells * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_nb,      (size_t)N * MAX_NEIGHBORS * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_nbCount,  N * sizeof(int)));

    // ---- Upload particle data (AoS → SoA) ----
    double *h_x  = new double[N], *h_y  = new double[N], *h_z  = new double[N];
    double *h_vx = new double[N], *h_vy = new double[N], *h_vz = new double[N];

    for (int i = 0; i < N; i++) {
        h_x[i]  = p[i].x;   h_y[i]  = p[i].y;   h_z[i]  = p[i].z;
        h_vx[i] = p[i].vx;  h_vy[i] = p[i].vy;  h_vz[i] = p[i].vz;
    }

    CUDA_CHECK(cudaMemcpy(d_x,  h_x,  N*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_y,  h_y,  N*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_z,  h_z,  N*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vx, h_vx, N*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vy, h_vy, N*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vz, h_vz, N*sizeof(double), cudaMemcpyHostToDevice));

    // ---- Initial neighbor list build + forces ----
    rebuildNeighborList(
        d_x, d_y, d_z,
        d_cellIDs, d_sortedIdx, d_cellStart, d_cellEnd,
        d_nb, d_nbCount, d_px, d_py, d_pz,
        N, nc, numCells, cellSize, L, rc_skin2, numBlocks);

    CUDA_CHECK(cudaMemset(d_fx, 0, N*sizeof(double)));
    CUDA_CHECK(cudaMemset(d_fy, 0, N*sizeof(double)));
    CUDA_CHECK(cudaMemset(d_fz, 0, N*sizeof(double)));
    CUDA_CHECK(cudaMemset(d_U,  0, N*sizeof(double)));

    ljForceKernel<<<numBlocks, GPU_BLOCK_SIZE>>>(
        d_x, d_y, d_z, d_fx, d_fy, d_fz, d_U,
        d_nb, d_nbCount, N, MAX_NEIGHBORS,
        sp->epsilon, sp->sigma, L, sp->rc2);
    CUDA_CHECK(cudaDeviceSynchronize());

    // ---- Main loop ----
    int rebuildCount = 0;
    auto t_start = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < nsteps; step++) {

        nbl_verletStage1<<<numBlocks, GPU_BLOCK_SIZE>>>(
            d_x, d_y, d_z, d_vx, d_vy, d_vz,
            d_fx, d_fy, d_fz, sp->dt, L, N);
        CUDA_CHECK(cudaDeviceSynchronize());

        if (needsRebuild(d_x, d_y, d_z, d_px, d_py, d_pz,
                         d_disp2, L, skin, N, numBlocks)) {
            rebuildNeighborList(
                d_x, d_y, d_z,
                d_cellIDs, d_sortedIdx, d_cellStart, d_cellEnd,
                d_nb, d_nbCount, d_px, d_py, d_pz,
                N, nc, numCells, cellSize, L, rc_skin2, numBlocks);
            rebuildCount++;
        }

        CUDA_CHECK(cudaMemset(d_fx, 0, N*sizeof(double)));
        CUDA_CHECK(cudaMemset(d_fy, 0, N*sizeof(double)));
        CUDA_CHECK(cudaMemset(d_fz, 0, N*sizeof(double)));
        CUDA_CHECK(cudaMemset(d_U,  0, N*sizeof(double)));

        ljForceKernel<<<numBlocks, GPU_BLOCK_SIZE>>>(
            d_x, d_y, d_z, d_fx, d_fy, d_fz, d_U,
            d_nb, d_nbCount, N, MAX_NEIGHBORS,
            sp->epsilon, sp->sigma, L, sp->rc2);
        CUDA_CHECK(cudaDeviceSynchronize());

        nbl_verletStage2<<<numBlocks, GPU_BLOCK_SIZE>>>(
            d_vx, d_vy, d_vz, d_fx, d_fy, d_fz, sp->dt, N);
        CUDA_CHECK(cudaDeviceSynchronize());

        if (f_energy) {
            computeKE<<<numBlocks, GPU_BLOCK_SIZE>>>(d_vx, d_vy, d_vz, d_K, N);

            double K = thrust::reduce(thrust::device_ptr<double>(d_K),
                                      thrust::device_ptr<double>(d_K + N),
                                      0.0, thrust::plus<double>());
            double U = thrust::reduce(thrust::device_ptr<double>(d_U),
                                      thrust::device_ptr<double>(d_U + N),
                                      0.0, thrust::plus<double>());
            fprintf(f_energy, "%d,%.10f,%.10f,%.10f\n", step, K, U, K + U);
        }
    }

    cudaDeviceSynchronize();
    auto t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t_stop - t_start;

    // ---- Copy results back (SoA → AoS) ----
    CUDA_CHECK(cudaMemcpy(h_x,  d_x,  N*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_y,  d_y,  N*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_z,  d_z,  N*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_vx, d_vx, N*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_vy, d_vy, N*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_vz, d_vz, N*sizeof(double), cudaMemcpyDeviceToHost));

    for (int i = 0; i < N; i++) {
        p[i].x  = h_x[i];   p[i].y  = h_y[i];   p[i].z  = h_z[i];
        p[i].vx = h_vx[i];  p[i].vy = h_vy[i];  p[i].vz = h_vz[i];
    }

    // ---- Cleanup ----
    cudaFree(d_x);  cudaFree(d_y);  cudaFree(d_z);
    cudaFree(d_vx); cudaFree(d_vy); cudaFree(d_vz);
    cudaFree(d_fx); cudaFree(d_fy); cudaFree(d_fz);
    cudaFree(d_U);  cudaFree(d_K);
    cudaFree(d_px); cudaFree(d_py); cudaFree(d_pz);
    cudaFree(d_disp2);
    cudaFree(d_cellIDs); cudaFree(d_sortedIdx);
    cudaFree(d_cellStart); cudaFree(d_cellEnd);
    cudaFree(d_nb); cudaFree(d_nbCount);

    delete[] h_x;  delete[] h_y;  delete[] h_z;
    delete[] h_vx; delete[] h_vy; delete[] h_vz;

    return elapsed.count();
}
