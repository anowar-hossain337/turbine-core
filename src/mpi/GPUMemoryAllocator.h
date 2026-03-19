
#pragma once

#ifndef TURBINECORE_GPUMEMORYALLOCATOR_H
#define TURBINECORE_GPUMEMORYALLOCATOR_H

#include "wind_turbine_core/ProjectDefines.h"

namespace turbine_core {

    namespace mpi {

        struct GPUHostMemoryAllocator
        {
            static void *allocate( size_t size )
            {
                void *p;
                WALBERLA_GPU_CHECK( cudaMallocHost( &p, size ))
                return p;
            }

            static void deallocate( void *ptr )
            {
                WALBERLA_GPU_CHECK( cudaFreeHost( ptr ))
            }

            static void memcpy( void *dst, void *src, size_t count )
            {
                std::memcpy( dst, src, count );
            }
        };

        struct GPUDeviceMemoryAllocator
        {
            static void *allocate( size_t size )
            {
                void *p;
                WALBERLA_GPU_CHECK( cudaMalloc( &p, size ));
                return p;
            }

            static void deallocate( void *ptr )
            {
                WALBERLA_GPU_CHECK( cudaFree( ptr ));
            }

            static void memcpy( void *dst, void *src, size_t count )
            {
                cudaMemcpy( dst, src, count, cudaMemcpyDeviceToDevice );
            }
        };

    }

}

#endif //TURBINECORE_GPUMEMORYALLOCATOR_H
