
#pragma once

#ifndef TURBINECORE_CONTROLMODEL_SimpleYawController_H
#define TURBINECORE_CONTROLMODEL_SimpleYawController_H

#include "ControlModel.h"
#include "conversion/Conversion.h"

namespace turbine_core {

    namespace control_model {

        class SimpleYawController final: public ControlModel {

            friend class ControlModel;

        public:

            static constexpr Type TYPE = Type::SIMPLE_YAW_CONTROLLER;

            HOST_DEVICE_PREFIX SimpleYawController() {}

            HOST_DEVICE_PREFIX SimpleYawController( const real_t & yawingSpeed)
            {
                yawingSpeed_ = yawingSpeed;
            }

            HOST_DEVICE_PREFIX SimpleYawController( const SimpleYawController & yawController)
                    : yawingSpeed_(yawController.yawingSpeed_) {}

            HOST_DEVICE_PREFIX SimpleYawController & operator=( SimpleYawController && model ) = delete;

            HOST_DEVICE_PREFIX ~SimpleYawController() override = default;

            HOST_DEVICE_PREFIX auto & yawingSpeed() {
                return yawingSpeed_;
            }

        private:

            HOST_DEVICE_PREFIX SimpleYawController * do_clone() const override {
                return new SimpleYawController(*this);
            }

            HOST_DEVICE_PREFIX void do_print() const override {
                printf("simpleYawController Control Model\n");
            }

            HOST_DEVICE_PREFIX void do_applyControl(Discretisation * discretisation, const uint_t timeStep, 
                                                    const real_t rotorMeanWind, const real_t rotorTorque, 
                                                    const real_t rotorOmega, const Vector3<real_t> rotorWindVaneWind) override {

                auto orientation = discretisation->points()[0].orientation;
                auto localVelocity = orientation.rotate(rotorWindVaneWind); // Velocity along the hub axis

                Vector3<real_t> projectedVelocity = orientation.getInverse().rotate(rotorWindVaneWind);
                real_t nacelleYaw = std::atan2(-projectedVelocity[1], -projectedVelocity[2]);

                if(nacelleYaw > 0) {
                    Vector3<real_t> rotSpeed(-yawingSpeed_, 0, 0);
                    discretisation->setRelativeRotationalVelocity( rotSpeed );
                } else {
                    Vector3<real_t> rotSpeed(+yawingSpeed_, 0, 0);
                    discretisation->setRelativeRotationalVelocity( rotSpeed );
                }
            }

            HOST_DEVICE_PREFIX void do_checkNeeds(bool & needsMeanVelocity, bool & needsTorque, bool & needsWindVane) override {
                needsMeanVelocity = false;
                needsTorque = false;
                needsWindVane = true;
                return;
            };

            HOST_DEVICE_PREFIX walberla::mpi::MPISize do_getItemSize() override {
                return sizeof(real_t);
            };

            HOST_DEVICE_PREFIX void do_fillBuffer(BufferType_T * buffer) override {

                memcpy(buffer, &yawingSpeed_, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                return;
            };

            HOST_DEVICE_PREFIX void do_unfillBuffer(BufferType_T * buffer) override {

                memcpy(&yawingSpeed_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                return;
            };

            real_t yawingSpeed_{0.};
        };

    } // namespace control_model

} // namespace turbine_core

#endif //TURBINECORE_CONTROLMODEL_RAWSAngles_H
