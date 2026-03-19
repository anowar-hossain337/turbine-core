
#pragma once

#ifndef TURBINECORE_GEOMETRICPOINT_H
#define TURBINECORE_GEOMETRICPOINT_H

#include <core/DataTypes.h>

#include "wind_turbine_core/math/Vector3.h"
#include "wind_turbine_core/math/Quaternion.h"

namespace turbine_core {

    namespace points {

        template<typename T>
        struct Point3 {

            Vector3 <T> position{};                 // global coordinates
            Vector3 <T> relativePosition{};         // parent-local coordinates
            Vector3 <T> velocity{};                 // global velocity
            Vector3 <T> rotationalVelocity{};       // global rotational velocity

            Quaternion<T> relativeOrientation{};    // rotation that transforms local to parent-local coordinates
            Quaternion<T> orientation{};            // rotation that transforms local to global coordinates

            HOST_DEVICE_PREFIX void updatePoint(const Point3 & reference, const Quaternion<T> & fullRotationMatrix,
                                                const Vector3 <T> & relativeRotationalVelocity = Vector3<T>()) {

                orientation = reference.orientation * relativeOrientation * fullRotationMatrix;

                auto globalRelativePosition = reference.orientation.rotate(relativePosition);

                position = reference.position + globalRelativePosition;
                velocity = reference.velocity + reference.rotationalVelocity % globalRelativePosition; // do not account for point.rotationalVelocity as distance is 0
                rotationalVelocity = reference.rotationalVelocity + orientation.rotate(relativeRotationalVelocity);

            }

        };

        //TODO assertion for valid size
        template<typename Buffer_T, typename Type_T>
        HOST_DEVICE_PREFIX Buffer_T * operator<<(Buffer_T * buffer, const Point3<Type_T> & point) {

            buffer << point.position;
            buffer << point.relativePosition;
            buffer << point.velocity;
            buffer << point.rotationalVelocity;

            buffer << point.relativeOrientation;
            buffer << point.orientation;

            return buffer;
        }

        //TODO assertion for valid size
        template<typename Buffer_T, typename Type_T>
        HOST_DEVICE_PREFIX Buffer_T * operator>>(Buffer_T * buffer, Point3<Type_T> & point) {

            buffer >> point.position;
            buffer >> point.relativePosition;
            buffer >> point.velocity;
            buffer >> point.rotationalVelocity;

            buffer >> point.relativeOrientation;
            buffer >> point.orientation;

            return buffer;
        }

    }

    using points::Point3;

}

#endif //TURBINECORE_GEOMETRICPOINT_H
