#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <chrono>

#include "gpu_manhattan_cell_cap.h"
#include "config.h"

// Must match CELL_CAP in cell_list.h
#ifndef CELL_CAP
#define CELL_CAP 16
#endif

static inline void ck(cudaError_t e, const char *msg) {
    if (e != cudaSuccess) {
        fprintf(stderr, "CUDA error (%s): %s\n", msg, cudaGetErrorString(e));
        std::exit(1);
    }
}

static int g_inited = 0;

// ---------------------------------------------------------------------------
// Device helpers
// ---------------------------------------------------------------------------

__device__ __forceinline__ int man_wrap(int a, int nc) {
    if (a < 0)   return a + nc;
    if (a >= nc) return a - nc;
    return a;
}

__device__ __forceinline__ int man_cell_index(int x, int y, int z, int nc) {
    return x + nc * (y + nc * z);
}

// ---------------------------------------------------------------------------
// Kernel: owner-computes-i, Manhattan (face-adjacent) neighbor cells only.
// Each thread accumulates force on particle i by scanning its own cell plus
// the 6 face-adjacent cells (±x, ±y, ±z).  Potential energy is double-counted
// here; the host multiplies by 0.5.
// ---------------------------------------------------------------------------
__global__ void manhattan_forces_U_kernel(
    const float * __restrict__ x,
    const float * __restrict__ y,
    const float * __restrict__ z,
    float       * __restrict__ fx,
    float       * __restrict__ fy,
    float       * __restrict__ fz,
    float       * __restrict__ Ui,
    const int   * __restrict__ counts,
    const int   * __restrict__ cells,
    unsigned N,
    int    nc,
    float  L,
    float  rc2,
    float  sigma,
    float  epsilon,
    float  cell_size
) {
    unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    float xi = x[i], yi = y[i], zi = z[i];

    int cx = man_wrap((int)floorf(xi / cell_size), nc);
    int cy = man_wrap((int)floorf(yi / cell_size), nc);
    int cz = man_wrap((int)floorf(zi / cell_size), nc);

    // 7 cells: self + 6 face-adjacent
    int neigh[7];
    neigh[0] = man_cell_index(cx,                  cy,                  cz,                  nc);
    neigh[1] = man_cell_index(man_wrap(cx + 1, nc), cy,                  cz,                  nc);
    neigh[2] = man_cell_index(man_wrap(cx - 1, nc), cy,                  cz,                  nc);
    neigh[3] = man_cell_index(cx,                  man_wrap(cy + 1, nc), cz,                  nc);
    neigh[4] = man_cell_index(cx,                  man_wrap(cy - 1, nc), cz,                  nc);
    neigh[5] = man_cell_index(cx,                  cy,                  man_wrap(cz + 1, nc), nc);
    neigh[6] = man_cell_index(cx,                  cy,                  man_wrap(cz - 1, nc), nc);

    float fxi = 0.f, fyi = 0.f, fzi = 0.f, ui = 0.f;

    for (int n = 0; n < 7; n++) {
        int c    = neigh[n];
        int cnt  = counts[c];
        int base = c * CELL_CAP;

        for (int k = 0; k < cnt; k++) {
            int j = cells[base + k];
            if (j < 0 || (unsigned)j == i) continue;

            float dx = xi - x[j];
            float dy = yi - y[j];
            float dz = zi - z[j];

            // minimum image
            dx -= L * nearbyintf(dx / L);
            dy -= L * nearbyintf(dy / L);
            dz -= L * nearbyintf(dz / L);

            float r2 = dx*dx + dy*dy + dz*dz;
            if (r2 > rc2 || r2 < 1e-20f) continue;

            float inv_r2  = 1.0f / r2;
            float sig2_r2 = (sigma * sigma) * inv_r2;
            float sr6     = sig2_r2 * sig2_r2 * sig2_r2;
            float sr12    = sr6 * sr6;

            ui += 4.0f * epsilon * (sr12 - sr6);

            float F_over_r = 24.0f * epsilon * (2.0f * sr12 - sr6) * inv_r2;
            fxi += F_over_r * dx;
            fyi += F_over_r * dy;
            fzi += F_over_r * dz;
        }
    }

    fx[i] = fxi;
    fy[i] = fyi;
    fz[i] = fzi;
    Ui[i] = ui;
}

// ---------------------------------------------------------------------------
// Low-level GPU interface
// ---------------------------------------------------------------------------

int gpu_manhattan_cuda_init(void) {
    if (g_inited) return 1;
    int ndev = 0;
    if (cudaGetDeviceCount(&ndev) != cudaSuccess || ndev <= 0) return 0;
    g_inited = 1;
    return 1;
}

void gpu_manhattan_cuda_cleanup(void) {
    g_inited = 0;
}

double gpu_manhattan_cuda_forces_and_U(Particle *p, const SimParams *sp, const CellList *cl) {
    if (!g_inited) return 0.0;

    const unsigned N = (unsigned)sp->N;
    const int ncell  = cl->ncell;

    std::vector<float> hx(N), hy(N), hz(N);
    for (unsigned i = 0; i < N; i++) {
        hx[i] = (float)p[i].x;
        hy[i] = (float)p[i].y;
        hz[i] = (float)p[i].z;
    }

    float *dx = nullptr, *dy = nullptr, *dz = nullptr;
    float *dfx = nullptr, *dfy = nullptr, *dfz = nullptr, *dUi = nullptr;
    int   *dcounts = nullptr, *dcells = nullptr;
    const size_t cells_bytes = (size_t)ncell * CELL_CAP * sizeof(int);

    ck(cudaMalloc(&dx,      sizeof(float) * N), "malloc x");
    ck(cudaMalloc(&dy,      sizeof(float) * N), "malloc y");
    ck(cudaMalloc(&dz,      sizeof(float) * N), "malloc z");
    ck(cudaMalloc(&dfx,     sizeof(float) * N), "malloc fx");
    ck(cudaMalloc(&dfy,     sizeof(float) * N), "malloc fy");
    ck(cudaMalloc(&dfz,     sizeof(float) * N), "malloc fz");
    ck(cudaMalloc(&dUi,     sizeof(float) * N), "malloc Ui");
    ck(cudaMalloc(&dcounts, sizeof(int) * ncell), "malloc counts");
    ck(cudaMalloc(&dcells,  cells_bytes),         "malloc cells");

    ck(cudaMemcpy(dx,      hx.data(),   sizeof(float) * N, cudaMemcpyHostToDevice), "cpy x");
    ck(cudaMemcpy(dy,      hy.data(),   sizeof(float) * N, cudaMemcpyHostToDevice), "cpy y");
    ck(cudaMemcpy(dz,      hz.data(),   sizeof(float) * N, cudaMemcpyHostToDevice), "cpy z");
    ck(cudaMemcpy(dcounts, cl->counts,  sizeof(int) * ncell, cudaMemcpyHostToDevice), "cpy counts");
    ck(cudaMemcpy(dcells,  cl->cells,   cells_bytes, cudaMemcpyHostToDevice),         "cpy cells");

    const int threads = GPU_BLOCK_SIZE;
    const int blocks  = ((int)N + threads - 1) / threads;

    manhattan_forces_U_kernel<<<blocks, threads>>>(
        dx, dy, dz, dfx, dfy, dfz, dUi,
        dcounts, dcells,
        N, cl->nc,
        (float)sp->L, (float)sp->rc2,
        (float)sp->sigma, (float)sp->epsilon,
        (float)cl->cell_size
    );
    ck(cudaGetLastError(),   "kernel launch");
    ck(cudaDeviceSynchronize(), "sync");

    std::vector<float> hfx(N), hfy(N), hfz(N), hUi(N);
    ck(cudaMemcpy(hfx.data(), dfx, sizeof(float) * N, cudaMemcpyDeviceToHost), "cpy fx");
    ck(cudaMemcpy(hfy.data(), dfy, sizeof(float) * N, cudaMemcpyDeviceToHost), "cpy fy");
    ck(cudaMemcpy(hfz.data(), dfz, sizeof(float) * N, cudaMemcpyDeviceToHost), "cpy fz");
    ck(cudaMemcpy(hUi.data(), dUi, sizeof(float) * N, cudaMemcpyDeviceToHost), "cpy Ui");

    for (unsigned i = 0; i < N; i++) {
        p[i].fx = (double)hfx[i];
        p[i].fy = (double)hfy[i];
        p[i].fz = (double)hfz[i];
    }

    // Ui is double-counted (each pair appears once from i's side and once from j's)
    double U = 0.0;
    for (unsigned i = 0; i < N; i++) U += (double)hUi[i];
    U *= 0.5;

    cudaFree(dx);  cudaFree(dy);  cudaFree(dz);
    cudaFree(dfx); cudaFree(dfy); cudaFree(dfz); cudaFree(dUi);
    cudaFree(dcounts); cudaFree(dcells);

    return U;
}

// ---------------------------------------------------------------------------
// High-level wrapper — same interface as cuda_full / cuda_cell_list / cuda_nbl
// Integration (velocity Verlet) runs on CPU; force evaluation runs on GPU.
// ---------------------------------------------------------------------------

static inline double wrap_pos(double x, double L) {
    x = fmod(x, L);
    if (x < 0) x += L;
    return x;
}

static double kinetic_energy_cpu(const Particle *p, size_t N) {
    double K = 0.0;
    for (size_t i = 0; i < N; i++)
        K += 0.5 * (p[i].vx*p[i].vx + p[i].vy*p[i].vy + p[i].vz*p[i].vz);
    return K;
}

double cuda_manhattan(const SimParams *sp, Particle *p, FILE *f_energy, int nsteps) {
    if (!gpu_manhattan_cuda_init()) {
        fprintf(stderr, "cuda_manhattan: CUDA init failed\n");
        return 0.0;
    }

    CellList cl;
    cell_list_init(&cl, sp->N, sp->L, sp->rc);
    cell_list_build(&cl, p, sp->L);

    // Compute initial forces
    gpu_manhattan_cuda_forces_and_U(p, sp, &cl);

    const double dt = sp->dt;
    const double L  = sp->L;
    const size_t N  = sp->N;

    auto t0 = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < nsteps; step++) {
        // Half-kick + drift
        for (size_t i = 0; i < N; i++) {
            p[i].vx += 0.5 * dt * p[i].fx;
            p[i].vy += 0.5 * dt * p[i].fy;
            p[i].vz += 0.5 * dt * p[i].fz;

            p[i].x = wrap_pos(p[i].x + dt * p[i].vx, L);
            p[i].y = wrap_pos(p[i].y + dt * p[i].vy, L);
            p[i].z = wrap_pos(p[i].z + dt * p[i].vz, L);
        }

        // Rebuild cell list and compute forces on GPU
        cell_list_build(&cl, p, L);
        double U = gpu_manhattan_cuda_forces_and_U(p, sp, &cl);

        // Second half-kick
        for (size_t i = 0; i < N; i++) {
            p[i].vx += 0.5 * dt * p[i].fx;
            p[i].vy += 0.5 * dt * p[i].fy;
            p[i].vz += 0.5 * dt * p[i].fz;
        }

        if (f_energy) {
            double K = kinetic_energy_cpu(p, N);
            fprintf(f_energy, "%d,%.10f,%.10f,%.10f\n", step, K, U, K + U);
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t1 - t0;

    cell_list_free(&cl);
    gpu_manhattan_cuda_cleanup();
    return elapsed.count();
}
