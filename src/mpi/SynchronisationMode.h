
#pragma once

#ifndef TURBINECORE_SYNCHRONISATIONMODE_H
#define TURBINECORE_SYNCHRONISATIONMODE_H

#include <core/mpi/MPIManager.h>

namespace turbine_core {

    namespace mpi {

        enum SynchronisationMode {
            SYNC_ALL,
            SYNC_MACRO,
            SYNC_FORCE,
            SYNC_GEOMETRY,
            SYNC_CONTROL
        };

        HOST_DEVICE_PREFIX walberla::mpi::MPISize itemSize( const SynchronisationMode mode ) {
            switch (mode) {
                case SYNC_ALL:
                    assert(0); // TODO NOT IMPLEMENTED
                case SYNC_MACRO:
                    return 7 * sizeof(real_t);
                case SYNC_FORCE:
                    return 3 * sizeof(real_t);
                case SYNC_GEOMETRY:
                    return 4 * 3 * sizeof(real_t) + 2 * 4 * sizeof(real_t);
                default:
                    printf("Not implemented.\n");
                    assert(0);
            }
        }

    }

}

#endif //TURBINECORE_SYNCHRONISATIONMODE_H
