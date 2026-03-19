
#pragma once

#ifndef TURBINECORE_CPUMEMORYALLOCATOR_H
#define TURBINECORE_CPUMEMORYALLOCATOR_H

#include "wind_turbine_core/ProjectDefines.h"

namespace turbine_core {

    namespace mpi {

        struct CPUMemoryAllocator
        {
            HOST_PREFIX static void *allocate( size_t size )
            {
                void *p = (void*) std::malloc(size);
                return p;
            }

            HOST_PREFIX static void deallocate( void *ptr )
            {
                free(ptr);
            }

            HOST_PREFIX static void memcpy( void *dst, void *src, size_t count )
            {
                std::memcpy( dst, src, count );
            }
        };

    }

}

#endif //TURBINECORE_CPUMEMORYALLOCATOR_H
