
#pragma once

#ifndef TURBINECORE_SMOOTHEDDIRACDELTAKERNEL_H
#define TURBINECORE_SMOOTHEDDIRACDELTAKERNEL_H

#include <cmath>

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/math/Vector3.h"

namespace turbine_core {

    namespace projectors {

        struct SmoothedDiracDeltaKernel {

            static constexpr uint_t xWidth{1};
            static constexpr uint_t yWidth{1};
            static constexpr uint_t zWidth{1};

            HOST_DEVICE_PREFIX static inline real_t getWeight(const Vector3 <real_t> &X, const Vector3 <real_t> &x,
                                                              const real_t dx = real_t(1), const real_t dy = real_t(1), const real_t dz = real_t(1) ) {
                return getWeight(X[0], X[1], X[2], x[0], x[1], x[2],
                                 dx, dy, dz);
            }

            HOST_DEVICE_PREFIX static inline real_t getWeight(const real_t X, const real_t Y, const real_t Z,
                                                              const real_t x, const real_t y, const real_t z,
                                                              const real_t dx = real_t(1), const real_t dy = real_t(1), const real_t dz = real_t(1)) {

                return smoothedDeltaFunction((X-x)/dx) *
                       smoothedDeltaFunction((Y-y)/dy) *
                       smoothedDeltaFunction((Z-z)/dz);

            }

        private:

            HOST_DEVICE_PREFIX static inline real_t smoothedDeltaFunction(const real_t r) {
                real_t rAbs = fabs(r);
                if( rAbs <= real_t(0.5) ) {
                    return (real_t(1) + sqrt( - real_t(3) * r * r + real_t(1) ) ) * (real_t(1) / real_t(3));
                } else if ( rAbs < real_t(1.5) ) {
                    return (real_t(5) - real_t(3) * rAbs - sqrt( - real_t(3) * ( real_t(1) - rAbs ) * ( real_t(1) - rAbs ) + real_t(1) ) ) * ( real_t(1) / real_t(6) );
                } else {
                    return real_t(0);
                }

            }

        };

    }

}

#endif //TURBINECORE_SMOOTHEDDIRACDELTAKERNEL_H
