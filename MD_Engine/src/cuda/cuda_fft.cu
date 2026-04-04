#include <complex.h>
#include <cstdio>
#include <cstdlib>
#include <cuComplex.h>
#include <cufft.h>
#include <math.h>
#include <time.h>
// Assertion to check for errors
#define CUDA_SAFE_CALL(ans)                                                    \
  {                                                                            \
    gpuAssert((ans), (char *)__FILE__, __LINE__);                              \
  }
inline void gpuAssert(cudaError_t code, char *file, int line,
                      bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "CUDA_SAFE_CALL: %s %s %d\n", cudaGetErrorString(code),
            file, line);
    if (abort)
      exit(code);
  }
}

#define FFT_SIZE (2 << 26)

uint64_t threads = 16;
uint64_t threadsPerBlock = threads * threads;
uint64_t blocks = (FFT_SIZE + threadsPerBlock - 1) / threadsPerBlock;
#define TESTNAME "cuda"

// #define SAVE

void write2disk(double *in, uint32_t n, int type, const char *filename) {
  FILE *file = fopen(filename, "w");
  if (file == NULL) {
    printf("ERROR: Opening %s\n", filename);
    return;
  }

  double real, imag;
  // if type == 0, double
  // if type == 1, complex
  int incr = (type == 0) ? 1 : 2;

  for (uint64_t i = 0; i < n; i += incr) {
    real = in[i];
    imag = (type != 0) ? in[i + 1] : 0.0;
    fprintf(file, "%.9f + %.9fi\n", real, imag);
  }

  fclose(file);
}

__global__ void bit_reverse(cuDoubleComplex *data, uint32_t n, int log2n) {
  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  if (idx >= n)
    return;

  uint32_t x = idx;
  x = (x >> 16) | (x << 16);
  x = ((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8);
  x = ((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4);
  x = ((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2);
  x = ((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1);

  x = x >> (32 - log2n);

  // Only swap if idx < rev to avoid double-swapping
  if (idx < x) {
    cuDoubleComplex temp = data[idx];
    data[idx] = data[x];
    data[x] = temp;
  }
}

__global__ void fft_shared(cuDoubleComplex *data, uint32_t n) {

  __shared__ cuDoubleComplex shared[16 * 16];

  int tid = threadIdx.x;
  int gid = tid + blockIdx.x * blockDim.x;

  shared[tid] = data[gid];

  __syncthreads();

  // Stages 1 to 8 can be done with shared memory
  for (int stage = 1; stage < 9; stage++) {
    int m = 1 << stage;  // e.g., 2, 4, 8, ..., 256
    int half_m = m >> 1; // e.g., 1, 2, 4, ..., 128
    int j = tid % m;     // local position within span

    if (j < half_m) {
      double angle = -2.0 * M_PI * j / m;
      double wr, wi;
      sincos(angle, &wi, &wr);
      cuDoubleComplex w = make_cuDoubleComplex(wr, wi);
      cuDoubleComplex t = cuCmul(w, shared[tid + half_m]);
      cuDoubleComplex u = shared[tid];

      shared[tid] = cuCadd(u, t);
      shared[tid + half_m] = cuCsub(u, t);
    }
    // Synchronize every stage
    __syncthreads();
  }

  // Write back to global memory
  data[gid] = shared[tid];
}
__global__ void fft_stage(cuDoubleComplex *data, uint32_t n, int m) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;

  int half_m = m >> 1;

  int group = idx / half_m;
  int j = idx % half_m;

  uint64_t i0 = group * m + j;
  uint64_t i1 = i0 + half_m;

  double angle = -2.0 * M_PI * j / m;
  double wr, wi;
  sincos(angle, &wi, &wr); // note: sincos gives sin(angle), cos(angle)

  cuDoubleComplex w = make_cuDoubleComplex(wr, wi);
  cuDoubleComplex t = cuCmul(w, data[i1]);
  cuDoubleComplex u = data[i0];

  data[i0] = cuCadd(u, t);
  data[i1] = cuCsub(u, t);
}

__global__ void fft_shared_twiddle(cuDoubleComplex *data, uint32_t n,
                                   const double *reals, const double *imags) {

  __shared__ cuDoubleComplex shared[16 * 16];

  int tid = threadIdx.x;
  int gid = tid + blockIdx.x * blockDim.x;

  shared[tid] = data[gid];

  __syncthreads();

  // Stages 1 to 8 can be done with shared memory
  for (int stage = 1; stage < 9; stage++) {
    int m = 1 << stage;  // e.g., 2, 4, 8, ..., 256
    int half_m = m >> 1; // e.g., 1, 2, 4, ..., 128
    int j = tid % m;     // local position within span

    if (j < half_m) {
      uint64_t twiddle_idx = j * (n / m);
      cuDoubleComplex w =
          make_cuDoubleComplex(reals[twiddle_idx], imags[twiddle_idx]);

      cuDoubleComplex t = cuCmul(w, shared[tid + half_m]);
      cuDoubleComplex u = shared[tid];

      shared[tid] = cuCadd(u, t);
      shared[tid + half_m] = cuCsub(u, t);
    }
    // Synchronize every stage
    __syncthreads();
  }

  // Write back to global memory
  data[gid] = shared[tid];
}
__global__ void fft_stage_twiddle(cuDoubleComplex *data, uint32_t n, int m,
                                  const double *reals, const double *imags) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;

  int half_m = m >> 1;
  int group = idx / half_m;
  int j = idx % half_m;

  uint64_t i0 = group * m + j;
  uint64_t i1 = i0 + half_m;

  if (i1 >= n)
    return;

  uint64_t twiddle_idx = j * (n / m);

  cuDoubleComplex w =
      make_cuDoubleComplex(reals[twiddle_idx], imags[twiddle_idx]);

  cuDoubleComplex t = cuCmul(w, data[i1]);
  cuDoubleComplex u = data[i0];

  data[i0] = cuCadd(u, t);
  data[i1] = cuCsub(u, t);
}

double complexDistance(const double *a, const double *b, int index) {
  double realA = a[2 * index];
  double imagA = a[2 * index + 1];
  double realB = b[2 * index];
  double imagB = b[2 * index + 1];

  double realDiff = realA - realB;
  double imagDiff = imagA - imagB;

  return sqrt(realDiff * realDiff + imagDiff * imagDiff);
}

int main() {

#ifdef SAVE
  char file_path[50];
#endif

  // GPU Timing variables
  cudaEvent_t start, stop;
  cudaEvent_t plan_start, plan_stop;

  float elapsed_gpu, elapsed_plan;

  // Arrays on GPU global memory
  double *d_data;
  double *d_result;
  cuDoubleComplex *d_mine;

  // Arrays on the host memory
  double *h_data;
  double *h_result;
  double *h_result1;

  // Allocate GPU memory
  uint64_t allocSize = (FFT_SIZE) * sizeof(double) * 2;
  CUDA_SAFE_CALL(cudaMalloc((void **)&d_data, allocSize));
  CUDA_SAFE_CALL(cudaMalloc((void **)&d_result, allocSize));

  CUDA_SAFE_CALL(cudaMalloc((void **)&d_mine, allocSize));

  // Allocate arrays on host memory
  h_data = (double *)malloc(allocSize);
  h_result = (double *)malloc(allocSize);

  h_result1 = (double *)malloc(allocSize);
  for (uint32_t i = 0; i < FFT_SIZE * 2; i += 2) {
    h_data[i] = sin(2 * M_PI * i / 2 / FFT_SIZE);
    h_data[i + 1] = 0.0;
  }

  double *h_twiddleReal, *h_twiddleImag;
  double *d_twiddleReal, *d_twiddleImag;

  uint64_t twiddle_count = FFT_SIZE / 2; // One per half-m segment

  h_twiddleReal = (double *)malloc(sizeof(double) * twiddle_count);
  h_twiddleImag = (double *)malloc(sizeof(double) * twiddle_count);

  for (uint64_t j = 0; j < twiddle_count; j++) {
    double angle = -2.0 * M_PI * j / FFT_SIZE;
    sincos(angle, &h_twiddleImag[j], &h_twiddleReal[j]);
  }

  CUDA_SAFE_CALL(
      cudaMalloc((void **)&d_twiddleReal, sizeof(double) * twiddle_count));
  CUDA_SAFE_CALL(
      cudaMalloc((void **)&d_twiddleImag, sizeof(double) * twiddle_count));

  CUDA_SAFE_CALL(cudaMemcpy(d_twiddleReal, h_twiddleReal,
                            sizeof(double) * twiddle_count,
                            cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(d_twiddleImag, h_twiddleImag,
                            sizeof(double) * twiddle_count,
                            cudaMemcpyHostToDevice));

#ifdef SAVE
  sprintf(file_path, "data/%s/input_%d", TESTNAME, FFT_SIZE);
  write2disk(h_data, FFT_SIZE * 2, 1, file_path);
#endif

  // Select GPU
  CUDA_SAFE_CALL(cudaSetDevice(0));
  // Copy data to GPU

  // Create the cuda events
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventCreate(&plan_start);
  cudaEventCreate(&plan_stop);

  /* cuFFT */
  CUDA_SAFE_CALL(cudaMemcpy(d_data, h_data, FFT_SIZE * sizeof(double) * 2,
                            cudaMemcpyHostToDevice));

  // Record event on the default stream
  cudaEventRecord(plan_start, 0);

  // Create a cuFFT plan
  cufftHandle plan;
  cufftCreate(&plan);
  cufftPlan1d(&plan, FFT_SIZE, CUFFT_Z2Z, 1); // 1D Real-to-Complex transform
  cudaEventRecord(plan_stop, 0);
  cudaEventSynchronize(plan_stop);

  cudaEventRecord(start, 0);
  // Execute cuFFT plan (GPU-based FFT)
  cufftExecZ2Z(plan, (cufftDoubleComplex *)d_data,
               (cufftDoubleComplex *)d_result, CUFFT_FORWARD);
  // Stop and destroy the timer
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  cudaEventElapsedTime(&elapsed_gpu, start, stop);
  cudaEventElapsedTime(&elapsed_plan, plan_start, plan_stop);
  printf("\ncuFFT Implementation:\n");
  printf("GPU time: %f (msec), Plan time: %f (msec)\n", elapsed_gpu,
         elapsed_plan);
  // Copy results back to host
  CUDA_SAFE_CALL(cudaMemcpy(h_result, d_result, FFT_SIZE * sizeof(double) * 2,
                            cudaMemcpyDeviceToHost));

#ifdef SAVE
  sprintf(file_path, "data/%s/output_%d_0", TESTNAME, FFT_SIZE);
  write2disk(h_result, FFT_SIZE * 2, 1, file_path);
#endif

  /* OUR FFT */
  CUDA_SAFE_CALL(cudaMemcpy(d_mine, h_data, FFT_SIZE * sizeof(double) * 2,
                            cudaMemcpyHostToDevice));

  cudaEventRecord(plan_start, 0);
  int log2n = log2(FFT_SIZE); // Compute on host side

  bit_reverse<<<blocks, threadsPerBlock>>>(d_mine, FFT_SIZE, log2n);
  cudaDeviceSynchronize();
  cudaEventRecord(plan_stop, 0);

  /* Version 1*/
  blocks = (FFT_SIZE / 2 + threadsPerBlock - 1) / threadsPerBlock;
  cudaEventRecord(start, 0);
  for (int s = 2; s <= FFT_SIZE; s *= 2) {
    fft_stage<<<blocks, threadsPerBlock>>>(d_mine, FFT_SIZE, s);
  }

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsed_gpu, start, stop);
  cudaEventElapsedTime(&elapsed_plan, plan_start, plan_stop);

  printf("\nV1 CUDA FFT Implementation\n");
  printf("GPU time: %f (msec), Bit Reversal time: %f (msec)\n", elapsed_gpu,
         elapsed_plan);

  // Copy results back to host
  CUDA_SAFE_CALL(cudaMemcpy(h_result1, d_mine, FFT_SIZE * sizeof(double) * 2,
                            cudaMemcpyDeviceToHost));
  {
    double max = 0.0;
    double mean = 0.0;
    for (uint32_t i = 0; i < FFT_SIZE; i++) {
      double tmp = complexDistance(h_result, h_result1, i);
      if (tmp > max) {
        max = tmp;
      }
      mean += tmp;
    }

    printf("\nMax Distance: %.8f, Average Distance: %.8f\n", max,
           mean / FFT_SIZE);
  }

  /* Version 2*/
  CUDA_SAFE_CALL(cudaMemcpy(d_mine, h_data, FFT_SIZE * sizeof(double) * 2,
                            cudaMemcpyHostToDevice));
  blocks = (FFT_SIZE + threadsPerBlock - 1) / threadsPerBlock;
  bit_reverse<<<blocks, threadsPerBlock>>>(d_mine, FFT_SIZE, log2n);

  cudaEventRecord(start, 0);
  fft_shared<<<blocks, threadsPerBlock>>>(d_mine, FFT_SIZE);

  blocks = (FFT_SIZE / 2 + threadsPerBlock - 1) / threadsPerBlock;
  for (int s = 512; s <= FFT_SIZE; s *= 2) {
    fft_stage<<<blocks, threadsPerBlock>>>(d_mine, FFT_SIZE, s);
  }

  // Stop and destroy the timer
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsed_gpu, start, stop);
  cudaEventElapsedTime(&elapsed_plan, plan_start, plan_stop);

  printf("\nV2 CUDA FFT Implementation\n");
  printf("GPU time: %f (msec), Bit Reversal time: %f (msec)\n", elapsed_gpu,
         elapsed_plan);

  // Copy results back to host
  CUDA_SAFE_CALL(cudaMemcpy(h_result1, d_mine, FFT_SIZE * sizeof(double) * 2,
                            cudaMemcpyDeviceToHost));
  {
    double max = 0.0;
    double mean = 0.0;
    for (uint32_t i = 0; i < FFT_SIZE; i++) {
      double tmp = complexDistance(h_result, h_result1, i);
      if (tmp > max) {
        max = tmp;
      }
      mean += tmp;
    }

    printf("\nMax Distance: %.8f, Average Distance: %.8f\n", max,
           mean / FFT_SIZE);
  }
#ifdef SAVE
  sprintf(file_path, "data/%s/output_%d_1", TESTNAME, FFT_SIZE);
  write2disk(h_result1, FFT_SIZE * 2, 1, file_path);
#endif

  /* Version 3*/
  CUDA_SAFE_CALL(cudaMemcpy(d_mine, h_data, FFT_SIZE * sizeof(double) * 2,
                            cudaMemcpyHostToDevice));
  blocks = (FFT_SIZE + threadsPerBlock - 1) / threadsPerBlock;
  bit_reverse<<<blocks, threadsPerBlock>>>(d_mine, FFT_SIZE, log2n);

  cudaEventRecord(start, 0);
  blocks = (FFT_SIZE / 2 + threadsPerBlock - 1) / threadsPerBlock;
  for (int s = 2; s <= FFT_SIZE; s *= 2) {
    fft_stage_twiddle<<<blocks, threadsPerBlock>>>(
        d_mine, FFT_SIZE, s, d_twiddleReal, d_twiddleImag);
  }

  // Stop and destroy the timer
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsed_gpu, start, stop);
  cudaEventElapsedTime(&elapsed_plan, plan_start, plan_stop);

  printf("\nV3 CUDA FFT Implementation\n");
  printf("GPU time: %f (msec), Bit Reversal time: %f (msec)\n", elapsed_gpu,
         elapsed_plan);

  // Copy results back to host
  CUDA_SAFE_CALL(cudaMemcpy(h_result1, d_mine, FFT_SIZE * sizeof(double) * 2,
                            cudaMemcpyDeviceToHost));
  {
    double max = 0.0;
    double mean = 0.0;
    for (uint32_t i = 0; i < FFT_SIZE; i++) {
      double tmp = complexDistance(h_result, h_result1, i);
      if (tmp > max) {
        max = tmp;
      }
      mean += tmp;
    }

    printf("\nMax Distance: %.8f, Average Distance: %.8f\n", max,
           mean / FFT_SIZE);
  }
#ifdef SAVE
  sprintf(file_path, "data/%s/output_%d_1", TESTNAME, FFT_SIZE);
  write2disk(h_result1, FFT_SIZE * 2, 1, file_path);
#endif

  /* Version 4*/
  CUDA_SAFE_CALL(cudaMemcpy(d_mine, h_data, FFT_SIZE * sizeof(double) * 2,
                            cudaMemcpyHostToDevice));
  blocks = (FFT_SIZE + threadsPerBlock - 1) / threadsPerBlock;
  bit_reverse<<<blocks, threadsPerBlock>>>(d_mine, FFT_SIZE, log2n);

  cudaEventRecord(start, 0);
  fft_shared_twiddle<<<blocks, threadsPerBlock>>>(d_mine, FFT_SIZE,
                                                  d_twiddleReal, d_twiddleImag);

  blocks = (FFT_SIZE / 2 + threadsPerBlock - 1) / threadsPerBlock;
  for (int s = 512; s <= FFT_SIZE; s *= 2) {
    fft_stage_twiddle<<<blocks, threadsPerBlock>>>(
        d_mine, FFT_SIZE, s, d_twiddleReal, d_twiddleImag);
  }

  // Stop and destroy the timer
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsed_gpu, start, stop);
  cudaEventElapsedTime(&elapsed_plan, plan_start, plan_stop);

  printf("\nV4 CUDA FFT Implementation\n");
  printf("GPU time: %f (msec), Bit Reversal time: %f (msec)\n", elapsed_gpu,
         elapsed_plan);

  // Copy results back to host
  CUDA_SAFE_CALL(cudaMemcpy(h_result1, d_mine, FFT_SIZE * sizeof(double) * 2,
                            cudaMemcpyDeviceToHost));
  {
    double max = 0.0;
    double mean = 0.0;
    for (uint32_t i = 0; i < FFT_SIZE; i++) {
      double tmp = complexDistance(h_result, h_result1, i);
      if (tmp > max) {
        max = tmp;
      }
      mean += tmp;
    }

    printf("\nMax Distance: %.8f, Average Distance: %.8f\n", max,
           mean / FFT_SIZE);
  }
#ifdef SAVE
  sprintf(file_path, "data/%s/output_%d_1", TESTNAME, FFT_SIZE);
  write2disk(h_result1, FFT_SIZE * 2, 1, file_path);
#endif

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  cudaEventDestroy(plan_start);
  cudaEventDestroy(plan_stop);
  cufftDestroy(plan);
  cudaFree(d_data);
  cudaFree(d_result);
  cudaFree(d_twiddleReal);
  cudaFree(d_twiddleImag);
  free(h_twiddleReal);
  free(h_twiddleImag);
  free(h_data);
  free(h_result);
  free(h_result1);
  return 0;
}