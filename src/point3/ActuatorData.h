
#pragma once

#ifndef TURBINECORE_ACTUATORDATA_H
#define TURBINECORE_ACTUATORDATA_H

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

#include "wind_turbine_core/math/Vector3.h"

#include "aerodynamics/AirfoilPolar.h"

namespace turbine_core {

    namespace actuator {

        struct ActuatorData {

            real_t angleOfAttack{};
            real_t width{};

            real_t          density{};
            Vector3<real_t> windVelocityGlobal{};
            Vector3<real_t> velocityLocal{};
            Vector3<real_t> forcesGlobal{};

            aerodynamics::AirfoilPolar polar{};

            HOST_DEVICE_PREFIX ActuatorData() {}

            HOST_DEVICE_PREFIX ActuatorData( const real_t _width, const aerodynamics::AirfoilPolar & _polar )
            : width(_width), polar(_polar)
            {}

        };

    }

    using actuator::ActuatorData;

}

#endif //TURBINECORE_ACTUATORDATA_H
