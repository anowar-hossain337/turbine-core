
#pragma once

#ifndef TURBINECORE_DISK_H
#define TURBINECORE_DISK_H

#include "Discretisation.h"

#include <core/DataTypes.h>

#include "point3/Point3.h"

#include "wind_turbine_core/ProjectDefines.h"

namespace turbine_core {

    namespace discretisation {

        class Disk final : public Discretisation {

        public:

            HOST_DEVICE_PREFIX explicit Disk(const Vector3<real_t> & position)
                    : radius_{0}
            {
                point_.position = position;
            }

            HOST_DEVICE_PREFIX explicit Disk(const real_t radius, const Vector3<real_t> & position)
                    : radius_{radius}
            {
                point_.position = position;
            }

            HOST_DEVICE_PREFIX Disk(const real_t & radius, const Vector3<real_t> & relPos,
                                    const Vector3<real_t> & relRot, const Quaternion<real_t> & relOrient )
                    : radius_(radius), relativeRotationalVelocity_(relRot),
                      rotationMatrix_(Quaternion<real_t>::makeQuaternionFromXYZAngles(relativeRotationalVelocity_)),
                      fullRotationMatrix_(rotationMatrix_.getInverse())
            {
                point_.position = relPos; //TODO needed only for base -> is there a better way?
                point_.relativePosition = relPos;
                point_.relativeOrientation = relOrient;
            }

            HOST_DEVICE_PREFIX Disk(const Disk & disk)
                    : radius_(disk.radius_), point_(disk.point_),
                      relativeRotationalVelocity_(disk.relativeRotationalVelocity_),
                      rotationMatrix_(disk.rotationMatrix_), fullRotationMatrix_(disk.fullRotationMatrix_)
            {}

            HOST_DEVICE_PREFIX Disk( Disk && disk ) = delete;

            HOST_DEVICE_PREFIX Disk & operator=(const Disk & disk) {

                if(this != &disk) {
                    radius_ = disk.radius_;
                    point_ = disk.point_;
                    relativeRotationalVelocity_ = disk.relativeRotationalVelocity_;
                    rotationMatrix_ = disk.rotationMatrix_;
                    fullRotationMatrix_ = disk.fullRotationMatrix_;
                }

                return *this;
            }

            HOST_DEVICE_PREFIX Disk & operator=( Disk && disk ) = delete;

            HOST_DEVICE_PREFIX ~Disk() override {}

            HOST_DEVICE_PREFIX points::Point3<real_t> * points() override {
                return &point_;
            }

            HOST_DEVICE_PREFIX uint_t nPoints() override {
                return 1;
            }

            HOST_DEVICE_PREFIX auto relativeRotationalVelocity() const {
                return relativeRotationalVelocity_;
            }

            HOST_DEVICE_PREFIX auto radius() const {
                return radius_;
            }

        private:

            HOST_DEVICE_PREFIX Disk * do_clone() const override {
                return new Disk(*this);
            }

            HOST_DEVICE_PREFIX void do_print() const override {
                printf("Disk (r=%f) position = [%f, %f, %f]\n", radius_, point_.position[0], point_.position[1], point_.position[2]);
            }

            HOST_DEVICE_PREFIX void do_update( const Point3<real_t> * reference, const uint_t timestep ) override {

                const uint_t startUp = 100;

                auto rotInc = relativeRotationalVelocity_;

                if(timestep < startUp) {
                    real_t factor = real_t(0.5) * (real_t(1) - cos(real_t(timestep) / real_t(startUp) * math::pi));
                    rotInc = rotInc * factor;
                }
                rotationMatrix_ = Quaternion<real_t>::makeQuaternionFromXYZAngles(rotInc);

                fullRotationMatrix_ = fullRotationMatrix_ * rotationMatrix_;
                point_.updatePoint(*reference, fullRotationMatrix_, rotInc);

            }

            HOST_DEVICE_PREFIX Point3<real_t> * do_getReferencePointForChild() override {
                return &point_;
            }

            HOST_DEVICE_PREFIX void do_setRelativeRotationalVelocity( const Vector3<real_t> & relativeRotationalVelocity) override {
                relativeRotationalVelocity_ = relativeRotationalVelocity;
            }

            HOST_DEVICE_PREFIX void do_setRelativeTranslationVelocity( const Vector3<real_t> & relativeTranslationVelocity) override {
                return;
            }

            HOST_DEVICE_PREFIX void do_updateRelativeOrientation( const Quaternion<real_t> & relativeOrientation) override {
                point_.relativeOrientation = point_.relativeOrientation*relativeOrientation;
                return;
            }

            HOST_DEVICE_PREFIX const Quaternion<real_t> & do_getRelativeOrientation() const override {
                return point_.relativeOrientation;
            }


            real_t radius_{0};
            Point3<real_t> point_{};

            Vector3<real_t> relativeRotationalVelocity_{};
            Quaternion<real_t> rotationMatrix_{};
            Quaternion<real_t> fullRotationMatrix_{};

        };

    }

}

#endif //TURBINECORE_DISK_H
