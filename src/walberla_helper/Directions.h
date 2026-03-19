
#pragma once

#ifndef TURBINECORE_DIRECTIONS_H
#define TURBINECORE_DIRECTIONS_H

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

#include "wind_turbine_core/math/Vector3.h"

#include <stencil/Directions.h>

namespace turbine_core {

    namespace stencil {

        namespace internal {

            constexpr int NDIR = walberla::stencil::NR_OF_DIRECTIONS;
            constexpr walberla::stencil::Direction INVALID = walberla::stencil::Direction::INVALID_DIR;

            /// The x component for each direction  \ingroup stencil
            MANAGED_PREFIX DEVICE_PREFIX int cx[NDIR]{};

            /// The y component for each direction \ingroup stencil
            MANAGED_PREFIX DEVICE_PREFIX int cy[NDIR]{};

            /// The z component for each direction \ingroup stencil
            MANAGED_PREFIX DEVICE_PREFIX int cz[NDIR]{};

            MANAGED_PREFIX DEVICE_PREFIX walberla::stencil::Direction dir[NDIR]{
                INVALID, INVALID, INVALID, INVALID, INVALID, INVALID, INVALID, INVALID, INVALID,
                INVALID, INVALID, INVALID, INVALID, INVALID, INVALID, INVALID, INVALID, INVALID,
                INVALID, INVALID, INVALID, INVALID, INVALID, INVALID, INVALID, INVALID, INVALID
            };

        } // namespace internal

        template<typename Stencil_T>
        HOST_PREFIX FORCEINLINE void initStencilData() {
            for( uint_t q = 0; q < internal::NDIR; ++q ) {
                internal::cx[q] = walberla::stencil::cx[q];
                internal::cy[q] = walberla::stencil::cy[q];
                internal::cz[q] = walberla::stencil::cz[q];
            }

            for( uint_t q = 0; q < Stencil_T::Q; ++q ) {
                internal::dir[q] = Stencil_T::dir[q];
            }
        }

        template<typename Stencil_T>
        HOST_DEVICE_PREFIX int cx( const uint_t idx ) {
            return internal::cx[internal::dir[idx]];
        }

        template<typename Stencil_T>
        HOST_DEVICE_PREFIX int cy( const uint_t idx ) {
            return internal::cy[internal::dir[idx]];
        }

        template<typename Stencil_T>
        HOST_DEVICE_PREFIX int cz( const uint_t idx ) {
            return internal::cz[internal::dir[idx]];
        }

        template<typename Stencil_T>
        HOST_DEVICE_PREFIX Vector3<int> c( const uint_t idx ) {
            auto dir = internal::dir[idx];
            return Vector3<int>( internal::cx[dir], internal::cy[dir], internal::cz[dir] );
        }

    } // namespace stencil

} // namespace turbine_core

#endif //TURBINECORE_DIRECTIONS_H
