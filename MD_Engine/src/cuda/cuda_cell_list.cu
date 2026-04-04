// nvcc -arch=sm_80 cuda_cell_list.cu -o cuda_cell_list
// nvcc -O3 -arch=sm_80 --use_fast_math src/cuda/cuda_cell_list.cu src/pdb_importer.c src/md_io.c -I include -o md_sim


#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <thrust/sequence.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <vector>
#include <cstdio>
#include <chrono>

#include "md.h"
#include "pdb_importer.h"
#include "cuda_algs.h"
#include "config.h"
#include "md_io.h"

// Using double for high precision as per your typedef
struct Particle_cuda {
    double3 pos;   // x, y, z
    double3 vel;   // vx, vy, vz
    double3 force; // fx, fy, fz
};

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
    
    // Cutoff (r=2.5, r^2=6.25)
    if (distSq < 1e-30 || distSq > rc2) return {0, 0, 0};

    double inv_r2 = 1.0 / distSq;
    double sig2_r2 = (sigma*sigma) * inv_r2;
    double sr6  = sig2_r2 * sig2_r2 * sig2_r2;
    double sr12 = sr6 * sr6;

    double forceMag = 24.0 * epsilon * (2.0*sr12 - sr6) * inv_r2;

    return {diff.x * forceMag, diff.y * forceMag, diff.z * forceMag};
}

// 1. INTEGRATION STAGE 1
__global__ void verletStage1(Particle_cuda* p, double dt, double simSize, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    // x = x + v*dt + 0.5*a*dt^2
    double invM = 1.0; //assume mass 1 
    p[i].pos.x += p[i].vel.x * dt + 0.5 * (p[i].force.x * invM) * dt * dt;
    p[i].pos.y += p[i].vel.y * dt + 0.5 * (p[i].force.y * invM) * dt * dt;
    p[i].pos.z += p[i].vel.z * dt + 0.5 * (p[i].force.z * invM) * dt * dt;

    // Half-step velocity
    p[i].vel.x += 0.5 * (p[i].force.x * invM) * dt;
    p[i].vel.y += 0.5 * (p[i].force.y * invM) * dt;
    p[i].vel.z += 0.5 * (p[i].force.z * invM) * dt;

    // PBC Wrap
    if(p[i].pos.x < 0) p[i].pos.x += simSize; if(p[i].pos.x >= simSize) p[i].pos.x -= simSize;
    if(p[i].pos.y < 0) p[i].pos.y += simSize; if(p[i].pos.y >= simSize) p[i].pos.y -= simSize;
    if(p[i].pos.z < 0) p[i].pos.z += simSize; if(p[i].pos.z >= simSize) p[i].pos.z -= simSize;
}

// 2. CELL ID CALCULATION
// why only x???
__global__ void calculateCellIDs(Particle_cuda* p, int* cellIDs, int3 gridDim, double cellSize, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    int cx = (int)(p[i].pos.x / cellSize);
    int cy = (int)(p[i].pos.y / cellSize);
    int cz = (int)(p[i].pos.z / cellSize);

    cellIDs[i] = (cz * gridDim.y * gridDim.x) + (cy * gridDim.x) + cx;
}

// 3. FORCE KERNEL (3D Neighbor Search)
__global__ void computeForces3D(Particle_cuda* p, int* cellStart, int* cellEnd, int* d_indices, int3 gridDim,
     double simSize, double cellSize, int n, double sigma, double epsilon, double rc2) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
//    double3 myPos = p[i].pos;
//    double3 totalF = {0, 0, 0};

    int real_idx = d_indices[i]; 
    double3 myPos = p[real_idx].pos;
    double3 totalF = {0, 0, 0};
    if (i >= n) return;

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

                int neighborIdx = (nz * gridDim.y * gridDim.x) + (ny * gridDim.x) + nx;
                int start = cellStart[neighborIdx];
                int end = cellEnd[neighborIdx];

                if (start == -1) continue;

//                for (int j = start; j < end; ++j) {
//                    if (i == j) continue;
//                    double3 f = computeLJForce(myPos, p[j].pos, simSize, sigma, epsilon, rc2);
//                    totalF.x += f.x; totalF.y += f.y; totalF.z += f.z;
//                }

                for (int j_ptr = start; j_ptr < end; ++j_ptr) {
                    if (i == j_ptr) continue;
                    
                    // Look up the neighbor's REAL memory location
                    int neighbor_real_idx = d_indices[j_ptr]; 
                    double3 f = computeLJForce(myPos, p[neighbor_real_idx].pos, simSize, sigma, epsilon, rc2);
                    totalF.x += f.x; totalF.y += f.y; totalF.z += f.z;
                }
    
            }
        }
    }

    p[real_idx].force = totalF; // Write back to the original memory spot
//    p[i].force = totalF;
}

// 4. INTEGRATION STAGE 2
__global__ void verletStage2(Particle_cuda* p, double dt, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double invM = 1.0; //assume mass 1
    p[i].vel.x += 0.5 * (p[i].force.x * invM) * dt;
    p[i].vel.y += 0.5 * (p[i].force.y * invM) * dt;
    p[i].vel.z += 0.5 * (p[i].force.z * invM) * dt;
}

// Boundary Finder (from earlier)
__global__ void findCellBounds(int* cellIDs, int* cellStart, int* cellEnd, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    int current = cellIDs[i];
    if (i == 0) {
        cellStart[current] = 0;
    } else {
        int prev = cellIDs[i-1];
        if (current != prev) {
            cellEnd[prev] = i;
            cellStart[current] = i;
        }
    }
    if (i == n - 1) {
        cellEnd[current] = n;
    }
}

void initializeLattice(Particle_cuda* h_p, int n, double simSize) {
    // Determine how many particles per side to fit n particles in a cube
    int pSide = (int)ceil(pow((double)n, 1.0/3.0));
    double spacing = simSize / (double)pSide;

    printf("Initializing %d particles in a %dx%dx%d lattice using double3 structure...\n", n, pSide, pSide, pSide);

    for (int i = 0; i < n; i++) {
        // get 3D coordinates
        int ix = i % pSide;
        int iy = (i / pSide) % pSide;
        int iz = i / (pSide * pSide);

        // Assign to .pos double3
        h_p[i].pos.x = ix * spacing + (spacing * 0.5);
        h_p[i].pos.y = iy * spacing + (spacing * 0.5);
        h_p[i].pos.z = iz * spacing + (spacing * 0.5);

        // Initialize .vel double3
        h_p[i].vel.x = 0.0;
        h_p[i].vel.y = 0.0;
        h_p[i].vel.z = 0.0;

        // Initialize .force double3
        h_p[i].force.x = 0.0;
        h_p[i].force.y = 0.0;
        h_p[i].force.z = 0.0;

    }
}

void reportFinalState(Particle_cuda* h_p, int n, double simSize, double cellSize) {
    double totalKE = 0;
    double totalPE = 0;
    
    for (int i = 0; i < n; i++) {
        // 1. Kinetic Energy: 0.5 * m * v^2
        double v2 = h_p[i].vel.x * h_p[i].vel.x + 
                    h_p[i].vel.y * h_p[i].vel.y + 
                    h_p[i].vel.z * h_p[i].vel.z;
        totalKE += 0.5 * 1 * v2;

        // 2. Potential Energy (Lennard-Jones)
        // Note: For a true baseline, we just check a few particles or 
        // do an O(N^2) check here since it's only once at the end.
        for (int j = i + 1; j < n; j++) {
            double dx = h_p[j].pos.x - h_p[i].pos.x;
            double dy = h_p[j].pos.y - h_p[i].pos.y;
            double dz = h_p[j].pos.z - h_p[i].pos.z;
            
            // Basic PBC for the report
            if (dx >  simSize*0.5) dx -= simSize; if (dx < -simSize*0.5) dx += simSize;
            if (dy >  simSize*0.5) dy -= simSize; if (dy < -simSize*0.5) dy += simSize;
            if (dz >  simSize*0.5) dz -= simSize; if (dz < -simSize*0.5) dz += simSize;

            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 6.25 && r2 > 0.001) {
                double r6inv = 1.0 / (r2 * r2 * r2);
                totalPE += 4.0 * (r6inv * r6inv - r6inv);
            }
        }
    }

    printf("\n--- Simulation Status ---\n");
    printf("Particles:    %d\n", n);
    printf("Kinetic E:    %.6f\n", totalKE);
    printf("Potential E:  %.6f\n", totalPE);
    printf("Total Energy: %.6f\n", totalKE + totalPE);
    printf("-------------------------------\n");
}

void cuda_cell_list(const SimParams *sp, Particle_cuda *p){

    const int numParticles = sp->N;
    const double simSize = sp->L; // Small box for high density
    const double cellSize = sp->rc; 
    const int cellDim = simSize/cellSize;
    const int3 gridDim = { cellDim, cellDim, cellDim }; 
    const int numCells = cellDim*cellDim*cellDim;

//    Particle* h_particles = (Particle*)malloc(numParticles * sizeof(Particle));
//    // TODO: Initialize h_particles with mass=1.0 and random positions/vel
//
//    initializeLattice(h_particles, numParticles, simSize);
    struct Particle_cuda * h_particles = p;

    struct Particle_cuda *d_particles;
//    SimParams * d_sp;
    int *d_cellStart, *d_cellEnd, *d_cellIDs;

    cudaMalloc(&d_particles, numParticles * sizeof(Particle_cuda));
    cudaMalloc(&d_cellIDs, numParticles * sizeof(int));
    cudaMalloc(&d_cellStart, numCells * sizeof(int));
    cudaMalloc(&d_cellEnd, numCells * sizeof(int));

    int *d_indices;
    cudaMalloc(&d_indices, numParticles * sizeof(int));
    thrust::sequence(thrust::device, d_indices, d_indices + numParticles);

    cudaMemcpy(d_particles, h_particles, numParticles * sizeof(struct Particle_cuda), cudaMemcpyHostToDevice);

    int blockSize = 256;
    int numBlocks = (numParticles + blockSize - 1) / blockSize;

    reportFinalState(h_particles, numParticles, simSize, cellSize);

    auto start = std::chrono::high_resolution_clock::now();
    for (int step = 0; step < AUTOTUNE_N_TIMESTEPS; step++) {
        verletStage1<<<numBlocks, blockSize>>>(d_particles, sp->dt, simSize, numParticles);
        
        cudaMemset(d_cellStart, -1, numCells * sizeof(int));
        calculateCellIDs<<<numBlocks, blockSize>>>(d_particles, d_cellIDs, gridDim, cellSize, numParticles);
        
//        thrust::device_ptr<int> d_keys(d_cellIDs);
//        thrust::device_ptr<struct Particle_cuda> d_values(d_particles);
//        thrust::sort_by_key(d_keys, d_keys + numParticles, d_values);
        thrust::device_ptr<int> d_keys(d_cellIDs);      // The Cell IDs (Keys)
        thrust::device_ptr<int> d_values(d_indices);   // The Original Indices (Values)
            
        // This only moves 4-byte ints! Fast on A100.
        thrust::sort_by_key(d_keys, d_keys + numParticles, d_values);

        
        cudaMemset(d_cellStart, -1, numCells * sizeof(int));
        cudaMemset(d_cellEnd, -1, numCells * sizeof(int));
        findCellBounds<<<numBlocks, blockSize>>>(d_cellIDs, d_cellStart, d_cellEnd, numParticles);
        
        computeForces3D<<<numBlocks, blockSize>>>(d_particles, d_cellStart, d_cellEnd, d_indices, gridDim,
             simSize, cellSize, numParticles, sp->sigma, sp->epsilon,sp->rc2);
        
        verletStage2<<<numBlocks, blockSize>>>(d_particles, sp->dt, numParticles);
    }

    cudaDeviceSynchronize();
    auto end = std::chrono::high_resolution_clock::now();

    cudaMemcpy(h_particles, d_particles, numParticles * sizeof(struct Particle_cuda), cudaMemcpyDeviceToHost);
    // 2. Report Results
    std::chrono::duration<double> elapsed = end - start;
    printf("A100 Execution Time for %d steps: %.4f seconds\n", AUTOTUNE_N_TIMESTEPS, elapsed.count());
    reportFinalState(h_particles, numParticles, simSize, cellSize);

}

int main(int argc, char **argv){
    if (argc < 2) {
        printf("Usage: %s input_path\n", argv[0]);
        return 1;
    }

    const char *input_path = argv[1];

    printf("\n================ INITIALIZING SIMULATION ================\n");
    
    Particle * p0;
    SimParams sp;
    int success;
    success = read_input(input_path, &p0, &sp);

    sp.sigma    = CONF_SIGMA;
    sp.epsilon  = CONF_EPSILON;
    sp.rc       = CONF_CUTOFF;
    sp.rc2      = sp.rc * sp.rc;
    sp.dt       = CONF_DELTA_T;
    sp.nthreads = N_THREADS;

    cuda_cell_list(&sp, reinterpret_cast<Particle_cuda*>(p0));
    io_write_pdb("output/final_positions.pdb", p0, sp.N);

    return 0;

}