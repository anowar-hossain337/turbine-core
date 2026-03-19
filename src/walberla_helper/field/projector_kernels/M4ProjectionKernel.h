
#pragma once

#ifndef TURBINECORE_M4PROJECTIONKERNEL_H
#define TURBINECORE_M4PROJECTIONKERNEL_H

#include <cmath>

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/math/Vector3.h"

namespace turbine_core {

    namespace projectors {

        struct M4ProjectionKernel {

            // The number of ghost cells must be equal to the nearest integer of 2.*rWidth
            // When using a single ghost layer, rWidth should take a value of 0.75 at maximum
            static constexpr real_t rWidth{0.75};
            static constexpr uint_t xWidth{static_cast<uint_t>(2.0*rWidth)};
            static constexpr uint_t yWidth{static_cast<uint_t>(2.0*rWidth)};
            static constexpr uint_t zWidth{static_cast<uint_t>(2.0*rWidth)};

            static constexpr real_t sigma{2./3.};

            HOST_DEVICE_PREFIX static inline real_t getWeight(const Vector3 <real_t> &X, const Vector3 <real_t> &x,
                                                              const real_t dx = real_t(1), const real_t dy = real_t(1), const real_t dz = real_t(1) ) {
                return getWeight(X[0], X[1], X[2], x[0], x[1], x[2],
                                 dx, dy, dz);
            }

            HOST_DEVICE_PREFIX static inline real_t getWeight(const real_t X, const real_t Y, const real_t Z,
                                                              const real_t x, const real_t y, const real_t z,
                                                              const real_t dx = real_t(1), const real_t dy = real_t(1), const real_t dz = real_t(1)) {

                return M4KernelFunction((X-x)/dx) *
                       M4KernelFunction((Y-y)/dy) *
                       M4KernelFunction((Z-z)/dz);

            }

        private:

            HOST_DEVICE_PREFIX static inline real_t M4KernelFunction(const real_t r) {
                real_t qAbs = fabs(r/rWidth);
                if( qAbs < real_t(1.0) ) {
                    return sigma/rWidth * (real_t(0.25)*(2.-qAbs)*(2.-qAbs)*(2.-qAbs)-(1.-qAbs)*(1.-qAbs)*(1.-qAbs));
                } else if ( qAbs < real_t(2.0) ) {
                    return sigma/rWidth * (real_t(0.25)*(2.-qAbs)*(2.-qAbs)*(2.-qAbs));
                } else {
                    return real_t(0);
                }

            }

        };

    }

}

#endif //TURBINECORE_M4PROJECTIONKERNEL_H
