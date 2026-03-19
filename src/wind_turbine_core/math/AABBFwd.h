
#pragma once

#ifndef TURBINECORE_AABBFWD_H
#define TURBINECORE_AABBFWD_H

#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace math {

        template<typename Type_T>
        class GenericAABB;

        typedef GenericAABB<real_t> AABB;

    } // namespace math

    using math::AABB;

}

#endif //TURBINECORE_AABBFWD_H
