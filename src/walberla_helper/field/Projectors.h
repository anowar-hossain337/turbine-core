
#pragma once

#ifndef TURBINECORE_PROJECTORS_H
#define TURBINECORE_PROJECTORS_H

#include "wind_turbine_core/ProjectDefines.h"

namespace turbine_core {

    namespace projectors {

        //TODO better way than void*... I have  bad gut feeling about this (C++ Guidelines do not directly prohibit this, as input and output parameters are non-void*)

        /** \struct Projectors
         *  \brief Container for interpolators and distributor needed in Actuator Line and Actuator Disk kernels.
         *
         *  This struct allows for passing interpolators and distributors of a priori unknown data type to force kernels.
         *  These functions cannot be templated as they are eventually virtual.
         *  Also, a common base class for interpolators, and distributors does not exist. So we cannot pass by pointer to base class.
         */
        struct Projectors {

            template<typename DensityInterpolator_T, typename VelocityInterpolator_T, typename ForceDistributor_T>
            HOST_DEVICE_PREFIX Projectors( DensityInterpolator_T * densityInterpolator, VelocityInterpolator_T velocityInterpolator, ForceDistributor_T * forceDistributor, const uint_t nGhostLayers )
            : densityInterpolator_(densityInterpolator), velocityInterpolator_(velocityInterpolator),
              forceDistributor_(forceDistributor), nGhostLayers_(nGhostLayers)
            {}

            template<typename Interpolator_T>
            HOST_DEVICE_PREFIX auto getDensityInterpolator() const {
                return static_cast<Interpolator_T*>(densityInterpolator_);
            }

            template<typename Interpolator_T>
            HOST_DEVICE_PREFIX auto getVelocityInterpolator() const {
                return static_cast<Interpolator_T*>(velocityInterpolator_);
            }

            template<typename Distributor_T>
            HOST_DEVICE_PREFIX auto getForceDistributor() const {
                return static_cast<Distributor_T*>(forceDistributor_);
            }

            HOST_DEVICE_PREFIX uint_t nGhostLayers() const {
                return nGhostLayers_;
            }

        private:

            void * densityInterpolator_;
            void * velocityInterpolator_;
            void * forceDistributor_;

            const uint_t nGhostLayers_{};

        };

    }

}

#endif //TURBINECORE_PROJECTORS_H
