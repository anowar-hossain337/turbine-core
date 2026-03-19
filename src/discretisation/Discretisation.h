
#pragma once

#ifndef TURBINECORE_DISCRETISATION_H
#define TURBINECORE_DISCRETISATION_H

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/math/Vector3.h"
#include "wind_turbine_core/math/Quaternion.h"

#include <cassert>

namespace turbine_core {

    namespace points {
        template <typename T>
        class Point3;
    }

    namespace gpu {
        template<typename Discretisation_T>
        struct CopyToGPU;

        template<typename Discretisation_T>
        HOST_PREFIX void copy( CopyToGPU<Discretisation_T>& cpyToGPU, const std::shared_ptr<Discretisation_T> & cpuPtr );
    }

    namespace discretisation {

        class Discretisation {

        public:

            HOST_DEVICE_PREFIX virtual ~Discretisation() {}

            HOST_DEVICE_PREFIX virtual points::Point3<real_t> * points() = 0;
            HOST_DEVICE_PREFIX virtual uint_t nPoints() = 0;

            HOST_DEVICE_PREFIX void print() const {
                do_print();
            }

            HOST_DEVICE_PREFIX Discretisation * clone() const {
                return do_clone();
            }

            HOST_DEVICE_PREFIX void update( const points::Point3<real_t> * reference, const uint_t timestep ) {
                do_update(reference, timestep);
            }

            HOST_DEVICE_PREFIX auto getReferencePointForChild() {
                return do_getReferencePointForChild();
            }

            HOST_DEVICE_PREFIX void setRelativeRotationalVelocity( Vector3<real_t> & relativeRotationalVelocity) {
                return do_setRelativeRotationalVelocity(relativeRotationalVelocity);
            }

            HOST_DEVICE_PREFIX void setRelativeTranslationVelocity( Vector3<real_t> & relativeTranslationVelocity) {
                return do_setRelativeTranslationVelocity(relativeTranslationVelocity);
            }

            HOST_DEVICE_PREFIX auto getRelativeOrientation() const {
                return do_getRelativeOrientation();
            }

            HOST_DEVICE_PREFIX void updateRelativeOrientation( const Quaternion<real_t> & relativeOrientation) {
                return do_updateRelativeOrientation(relativeOrientation);
            }

        private:

            HOST_DEVICE_PREFIX virtual void do_print() const = 0;
            HOST_DEVICE_PREFIX virtual Discretisation * do_clone() const = 0;
            HOST_DEVICE_PREFIX virtual void do_update( const points::Point3<real_t> * reference, const uint_t timestep ) = 0;
            HOST_DEVICE_PREFIX virtual points::Point3<real_t> * do_getReferencePointForChild() = 0;
            HOST_DEVICE_PREFIX virtual void do_setRelativeRotationalVelocity(const Vector3<real_t> & relativeRotationalVelocity) = 0;
            HOST_DEVICE_PREFIX virtual void do_setRelativeTranslationVelocity(const Vector3<real_t> & relativeTranslationVelocity) = 0;
            HOST_DEVICE_PREFIX virtual const Quaternion<real_t> & do_getRelativeOrientation() const = 0;
            HOST_DEVICE_PREFIX virtual void do_updateRelativeOrientation( const Quaternion<real_t> & relativeOrientation) = 0;
        };

    }

    using discretisation::Discretisation;

}

#endif //TURBINECORE_DISCRETISATION_H
