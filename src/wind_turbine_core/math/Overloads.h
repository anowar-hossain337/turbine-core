
#pragma once

#ifndef TURBINECORE_OVERLOADS_H
#define TURBINECORE_OVERLOADS_H

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/ProjectDefines.h"

#include <cmath>

namespace turbine_core {

    namespace math {

        HOST_DEVICE_PREFIX FORCEINLINE
        real_t max(const real_t x, const real_t y) {
            return fmax(x, y);
        }

        HOST_DEVICE_PREFIX FORCEINLINE
        real_t min(const real_t x, const real_t y) {
            return fmin(x, y);
        }

        HOST_DEVICE_PREFIX FORCEINLINE
        uint_t max(const uint_t x, const uint_t y) {
            return x > y ? x : y;
        }

        HOST_DEVICE_PREFIX FORCEINLINE
        uint_t min(const uint_t x, const uint_t y) {
            return x > y ? y : x;
        }

    }
}


#endif //TURBINECORE_OVERLOADS_H
