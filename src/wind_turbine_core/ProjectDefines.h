
#pragma once

#ifndef TURBINECORE_PROJECTDEFINES_H
#define TURBINECORE_PROJECTDEFINES_H

#include <waLBerlaDefinitions.h>
#include <cstdio>

#ifdef WALBERLA_BUILD_WITH_CUDA

#include <cuda.h>
#include <cuda_runtime_api.h>

#define HOST_DEVICE_PREFIX __host__ __device__
#define HOST_PREFIX __host__
#define DEVICE_PREFIX __device__
#define GLOBAL_PREFIX __global__
#define CONSTANT_PREFIX __constant__
#define MANAGED_PREFIX __managed__
#define FORCEINLINE __forceinline__

#define LAUNCH_BOUNDS(maxThreadsPerBlock) __launch_bounds__(maxThreadsPerBlock)

#define ALIGN(n) __align__(n)

#ifndef NDEBUG
#define TURBINE_GPU_CHECK() { WALBERLA_GPU_CHECK(cudaPeekAtLastError()); WALBERLA_GPU_CHECK(cudaDeviceSynchronize()) }

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}
#else
#define TURBINE_GPU_CHECK()
#define gpuErrchk(ans)
#endif

#else

#define HOST_DEVICE_PREFIX
#define HOST_PREFIX
#define DEVICE_PREFIX
#define GLOBAL_PREFIX
#define CONSTANT_PREFIX
#define MANAGED_PREFIX
#define FORCEINLINE __inline__ __attribute__((always_inline))

#define LAUNCH_BOUNDS(maxThreadsPerBlock)

#define ALIGN(n) __attribute__((aligned(n)))

#endif

#define RESTRICT __restrict__

#endif //TURBINECORE_PROJECTDEFINES_H
