#ifndef COMMON_H_
#define COMMON_H_

#include <sys/time.h>
#include <stdio.h>

#define CHECK(call) { \
  const cudaError_t error = call; \
  if (error != cudaSuccess) { \
    fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__); \
    fprintf(stderr, "code: %d, reason: %s\n", error, cudaGetErrorString(error)); \
  } \
}

#define CHECK_CUBLAS(call) { \
  cublasStatus_t err; \
  if ((err = (call)) != CUBLAS_STATUS_SUCCESS) { \
    fprintf(stderr, "Got CUBLAS error %d at %s:%d\n", err, __FILE__, __LINE__); \
    exit(1); \
  } \
}

#define CHECK_CURAND(call) { \
  curandStatus_t err; \
  if ((err = (call)) != CURAND_STATUS_SUCCESS) { \
    fprintf(stderr, "Got CURAND error %d at %s:%d\n", err, __FILE__, __LINE__); \
    exit(1); \
  } \
}

#define CHECK_CUFFT(call) { \
  cufftResult err; \
  if ( (err = (call)) != CUFFT_SUCCESS) { \
    fprintf(stderr, "Got CUFFT error %d at %s:%d\n", err, __FILE__, __LINE__); \
    exit(1); \
  } \
}

#define CHECK_CUSPARSE(call) { \
  cusparseStatus_t err; \
  if ((err = (call)) != CUSPARSE_STATUS_SUCCESS) { \
    fprintf(stderr, "Got error %d at %s:%d\n", err, __FILE__, __LINE__); \
    cudaError_t cuda_err = cudaGetLastError(); \
    if (cuda_err != cudaSuccess) { \
      fprintf(stderr, "  CUDA error \"%s\" also detected\n", \
      cudaGetErrorString(cuda_err)); \
    } \
    exit(1); \
  } \
}

inline double seconds() {
  struct timeval tp;
  struct timezone tzp;
  int i = gettimeofday(&tp, &tzp);
  return (double)tp.tv_sec + (double)tp.tv_usec*1.e-6;
}

/*
 * Display a variety of information on the first CUDA device in this system,
 * including driver version, runtime version, compute capability, bytes of
 * global memory, etc.
 */

void printDevInfo(int dev) {
  int driverVersion = 0;
  int runtimeVersion = 0;
  cudaDeviceProp deviceProp;
  CHECK(cudaGetDeviceProperties(&deviceProp, dev));
  printf("Device %d: \"%s\"\n", dev, deviceProp.name);
  CHECK(cudaDriverGetVersion(&driverVersion));
  CHECK(cudaRuntimeGetVersion(&runtimeVersion));
  printf(
    "  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n",
    driverVersion / 1000, (driverVersion % 100) / 10,
    runtimeVersion / 1000, (runtimeVersion % 100) / 10);
  printf(
    "  Number of multiprocessors on device            %d\n",
    deviceProp.multiProcessorCount);
  printf(
    "  CUDA Capability Major/Minor version number:    %d.%d\n",
    deviceProp.major, deviceProp.minor);
  printf(
    "  Total amount of global memory:                 %.2f GBytes (%llu "
    "bytes)\n", (float)deviceProp.totalGlobalMem / pow(1024.0, 3),
    (unsigned long long)deviceProp.totalGlobalMem);
  printf(
    "  GPU Clock rate:                                %.0f MHz (%0.2f "
    "GHz)\n", deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);
  printf(
    "  Memory Clock rate:                             %.0f Mhz\n",
    deviceProp.memoryClockRate * 1e-3f);
  printf(
    "  Memory Bus Width:                              %d-bit\n",
    deviceProp.memoryBusWidth);
  if (deviceProp.l2CacheSize) {
    printf("  L2 Cache Size:                                 %d bytes\n",
    deviceProp.l2CacheSize);
  }
  printf(
    "  Max Texture Dimension Size (x,y,z)             1D=(%d), "
    "2D=(%d,%d), 3D=(%d,%d,%d)\n", deviceProp.maxTexture1D,
    deviceProp.maxTexture2D[0], deviceProp.maxTexture2D[1],
    deviceProp.maxTexture3D[0], deviceProp.maxTexture3D[1],
    deviceProp.maxTexture3D[2]);
  printf(
    "  Max Layered Texture Size (dim) x layers        1D=(%d) x %d, "
    "2D=(%d,%d) x %d\n", deviceProp.maxTexture1DLayered[0],
    deviceProp.maxTexture1DLayered[1], deviceProp.maxTexture2DLayered[0],
    deviceProp.maxTexture2DLayered[1], deviceProp.maxTexture2DLayered[2]);
  printf(
    "  Total amount of constant memory:               %lu bytes\n",
    deviceProp.totalConstMem);
  printf(
    "  Total amount of shared memory per block:       %lu bytes\n",
    deviceProp.sharedMemPerBlock);
  printf(
    "  Total number of registers available per block: %d\n",
    deviceProp.regsPerBlock);
  printf(
    "  Warp size:                                     %d\n",
    deviceProp.warpSize);
  printf(
    "  Maximum number of threads per multiprocessor:  %d\n",
    deviceProp.maxThreadsPerMultiProcessor);
  printf(
    "  Maximum number of threads per block:           %d\n",
    deviceProp.maxThreadsPerBlock);
  printf(
    "  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
    deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1],
    deviceProp.maxThreadsDim[2]);
  printf(
    "  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
    deviceProp.maxGridSize[0], deviceProp.maxGridSize[1],
    deviceProp.maxGridSize[2]);
  printf(
    "  Maximum memory pitch:                          %lu bytes\n",
    deviceProp.memPitch);
  printf("\n\n");
}
void printAllDevInfo() {
  int deviceCount;
  CHECK(cudaGetDeviceCount(&deviceCount));
  if (deviceCount == 0) {
    printf("There are no available device(s) that support CUDA\n");
  } else {
    printf("Detected %d CUDA Capable device(s)\n", deviceCount);
  }
  for (int dev = 0; dev < deviceCount; dev++) {
    printDevInfo(dev);
  }
}

int getOptDevByMulProcessorCount() {
  int numDevices = 0;
  CHECK(cudaGetDeviceCount(&numDevices));
  if (numDevices == 0) {
    printf("There are no available device(s) that support CUDA\n");
    exit(1);
  }
  printf("Detected %d CUDA Capable device(s)\n", numDevices);
  int maxMultiprocessors = 0;
  int maxDevice = 0;
  for (int dev = 0; dev < numDevices; dev++) {
    cudaDeviceProp props;
    CHECK(cudaGetDeviceProperties(&props, dev));
    if (maxMultiprocessors < props.multiProcessorCount) {
      maxMultiprocessors = props.multiProcessorCount;
      maxDevice = dev;
    }
  }
  printf("choose dev: %d\n", maxDevice);
  return maxDevice;
}

#endif // COMMON_H_
