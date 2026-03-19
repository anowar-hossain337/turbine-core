
#pragma once
#ifndef TURBINECORE_WALBERLADATATYPES_H
#define TURBINECORE_WALBERLADATATYPES_H

#include "wind_turbine_core/ProjectDefines.h"
#include <core/DataTypes.h>
#include <domain_decomposition/BlockDataID.h>

namespace turbine_core {

    using walberla::int32_t;
    using walberla::uint_t;
    using walberla::real_t;

    using walberla::cell_idx_t;

    using walberla::BlockDataID;

    namespace math {
        MANAGED_PREFIX DEVICE_PREFIX real_t pi{3.1415926535897932};
    }
}

#endif //TURBINECORE_WALBERLADATATYPES_H
