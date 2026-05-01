#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <chrono>
#include <cstdio>

#include "md.h"
#include "md_io.h"
#include "cuda_full.h"
#include "config.h"
#include "cuda_common.cuh"

// ---------------------------------------------------------------------------
// Kernels
// ---------------------------------------------------------------------------

// static to avoid link conflicts with nbl_verletStage1/2 under -rdc=true
static __global__ void fullVerletStage1(
    double *x,  double *y,  double *z,
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

    x[i] = fmod(x[i], L); if (x[i] < 0.0) x[i] += L;
    y[i] = fmod(y[i], L); if (y[i] < 0.0) y[i] += L;
    z[i] = fmod(z[i], L); if (z[i] < 0.0) z[i] += L;
}

static __global__ void fullVerletStage2(
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

static __global__ void fullComputeKE(
    const double *vx, const double *vy, const double *vz,
    double *K_per_particle, int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    K_per_particle[i] = 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
}

// O(N^2) all-pairs force with shared memory tiling.
// Full shell
// Each tile loads GPU_BLOCK_SIZE j-positions into shared memory so the
// whole block reuses them before fetching the next tile.
static __global__ void computeForcesFull(
    const double *x, const double *y, const double *z,
    double *fx, double *fy, double *fz,
    double *U_per_particle,
    int N, double L, double sigma, double epsilon, double rc2)
{
    __shared__ double sx[GPU_BLOCK_SIZE];
    __shared__ double sy[GPU_BLOCK_SIZE];
    __shared__ double sz[GPU_BLOCK_SIZE];

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    double xi = (i < N) ? x[i] : 0.0;
    double yi = (i < N) ? y[i] : 0.0;
    double zi = (i < N) ? z[i] : 0.0;

    double fx_i = 0.0, fy_i = 0.0, fz_i = 0.0, U_i = 0.0;
    double halfL = 0.5 * L;
    double sig2  = sigma * sigma;

    int numTiles = (N + GPU_BLOCK_SIZE - 1) / GPU_BLOCK_SIZE;

    for (int tile = 0; tile < numTiles; tile++) {
        int j_global = tile * GPU_BLOCK_SIZE + threadIdx.x;

        sx[threadIdx.x] = (j_global < N) ? x[j_global] : 0.0;
        sy[threadIdx.x] = (j_global < N) ? y[j_global] : 0.0;
        sz[threadIdx.x] = (j_global < N) ? z[j_global] : 0.0;
        __syncthreads();

        if (i < N) {
            for (int k = 0; k < GPU_BLOCK_SIZE; k++) {
                int j = tile * GPU_BLOCK_SIZE + k;
                if (j >= N || j == i) continue;

                double dx = xi - sx[k];
                double dy = yi - sy[k];
                double dz = zi - sz[k];

                if (dx >  halfL) dx -= L; else if (dx < -halfL) dx += L;
                if (dy >  halfL) dy -= L; else if (dy < -halfL) dy += L;
                if (dz >  halfL) dz -= L; else if (dz < -halfL) dz += L;

                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < 1e-30 || r2 > rc2) continue;

                double inv_r2  = 1.0 / r2;
                double sig2_r2 = sig2 * inv_r2;
                double sr6     = sig2_r2 * sig2_r2 * sig2_r2;
                double sr12    = sr6 * sr6;
                double fmag    = 24.0 * epsilon * (2.0*sr12 - sr6) * inv_r2;

                fx_i += fmag * dx;
                fy_i += fmag * dy;
                fz_i += fmag * dz;
                U_i  += 4.0 * epsilon * (sr12 - sr6);
            }
        }
        __syncthreads();
    }

    if (i < N) {
        fx[i]             = fx_i;
        fy[i]             = fy_i;
        fz[i]             = fz_i;
        U_per_particle[i] = 0.5 * U_i;
    }
}

// ---------------------------------------------------------------------------
//  entry point
// ---------------------------------------------------------------------------

double cuda_full(const SimParams *sp, Particle *p, FILE *f_energy, int nsteps)
{
    const int    N = (int)sp->N;
    const double L = sp->L;

    const int numBlocks = (N + GPU_BLOCK_SIZE - 1) / GPU_BLOCK_SIZE;

    // ---- Allocate device arrays (SoA) ----
    double *d_x,  *d_y,  *d_z;
    double *d_vx, *d_vy, *d_vz;
    double *d_fx, *d_fy, *d_fz;
    double *d_U,  *d_K;

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

    // ---- Upload (AoS -> SoA) ----
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

    // ---- Initial forces ----
    computeForcesFull<<<numBlocks, GPU_BLOCK_SIZE>>>(
        d_x, d_y, d_z, d_fx, d_fy, d_fz, d_U,
        N, L, sp->sigma, sp->epsilon, sp->rc2);
    CUDA_CHECK(cudaDeviceSynchronize());

    // ---- Main loop ----
    auto t_start = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < nsteps; step++) {
        fullVerletStage1<<<numBlocks, GPU_BLOCK_SIZE>>>(
            d_x, d_y, d_z, d_vx, d_vy, d_vz,
            d_fx, d_fy, d_fz, sp->dt, L, N);
        CUDA_CHECK(cudaDeviceSynchronize());

        computeForcesFull<<<numBlocks, GPU_BLOCK_SIZE>>>(
            d_x, d_y, d_z, d_fx, d_fy, d_fz, d_U,
            N, L, sp->sigma, sp->epsilon, sp->rc2);
        CUDA_CHECK(cudaDeviceSynchronize());

        fullVerletStage2<<<numBlocks, GPU_BLOCK_SIZE>>>(
            d_vx, d_vy, d_vz, d_fx, d_fy, d_fz, sp->dt, N);
        CUDA_CHECK(cudaDeviceSynchronize());

        if (f_energy) {
            fullComputeKE<<<numBlocks, GPU_BLOCK_SIZE>>>(d_vx, d_vy, d_vz, d_K, N);
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

    // ---- Copy results back (SoA -> AoS) ----
    CUDA_CHECK(cudaMemcpy(h_x,  d_x,  N*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_y,  d_y,  N*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_z,  d_z,  N*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_vx, d_vx, N*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_vy, d_vy, N*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_vz, d_vz, N*sizeof(double), cudaMemcpyDeviceToHost));

    for (int i = 0; i < N; i++) {
        p[i].x  = h_x[i];  p[i].y  = h_y[i];  p[i].z  = h_z[i];
        p[i].vx = h_vx[i]; p[i].vy = h_vy[i]; p[i].vz = h_vz[i];
    }

    // ---- Cleanup ----
    cudaFree(d_x);  cudaFree(d_y);  cudaFree(d_z);
    cudaFree(d_vx); cudaFree(d_vy); cudaFree(d_vz);
    cudaFree(d_fx); cudaFree(d_fy); cudaFree(d_fz);
    cudaFree(d_U);  cudaFree(d_K);

    delete[] h_x;  delete[] h_y;  delete[] h_z;
    delete[] h_vx; delete[] h_vy; delete[] h_vz;

    return elapsed.count();
}
