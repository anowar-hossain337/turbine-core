
#pragma once

#ifndef TURBINECORE_CONTROLMODEL_H
#define TURBINECORE_CONTROLMODEL_H

#include "mpi/CustomMemoryBuffer.h"
#include "mpi/CPUMemoryAllocator.h"

namespace turbine_core {

    namespace control_model {

        class ControlModel {

        protected:
            using BufferType_T = unsigned char;

        public:

            enum Type {
                NONE,
                RAWS_ANGLES,
                RAWS_OMEGA,
                TIME_CONTROL_VELOCITIES,
                DISCON_TORQUE_CONTROL,
                DISCON_PITCH_CONTROL,
                SIMPLE_YAW_CONTROLLER
            };

            static Type toType(const std::string & identifier) {

                auto id = string::toLowercase(identifier);
                id = string::removeSpaces(id);

                if (id == "none") {
                    return Type::NONE;
                } else if (id == "raws_angles") {
                    return Type::RAWS_ANGLES;
                } else if (id == "raws_omega") {
                    return Type::RAWS_OMEGA;
                } else if (id == "time_control_velocities") {
                    return Type::TIME_CONTROL_VELOCITIES;
                } else if (id == "discon_torque_control") {
                    return Type::DISCON_TORQUE_CONTROL;
                } else if (id == "discon_pitch_control") {
                    return Type::DISCON_PITCH_CONTROL;
                } else if (id == "simple_yaw_controller") {
                    return Type::SIMPLE_YAW_CONTROLLER;
                } else {
                    WALBERLA_ABORT("Invalid controller name (" << identifier << ").")
                }

            }

            static std::string toString(const ControlModel::Type & type) {

                switch (type) {
                    case Type::NONE :
                        return "None";
                    case Type::RAWS_ANGLES :
                        return "RAWS_Angles";
                    case Type::RAWS_OMEGA :
                        return "RAWS_Omega";
                    case Type::TIME_CONTROL_VELOCITIES :
                        return "timeControlVelocities";
                    case Type::DISCON_TORQUE_CONTROL :
                        return "DISCON_Torque_Control";
                    case Type::DISCON_PITCH_CONTROL :
                        return "DISCON_Pitch_Control";
                    case Type::SIMPLE_YAW_CONTROLLER :
                        return "simpleYawController";
                }

            }

            HOST_DEVICE_PREFIX ControlModel() {}

            HOST_DEVICE_PREFIX virtual ~ControlModel() {}

            HOST_DEVICE_PREFIX void print() const {
                do_print();
            }

            HOST_DEVICE_PREFIX ControlModel * clone() const {
                return do_clone();
            }

            HOST_DEVICE_PREFIX void applyControl(Discretisation * discretisation, const uint_t timeStep, const real_t rotorMeanWind, const real_t rotorTorque,
                const real_t rotorOmega, const Vector3<real_t> rotorWindVaneWind) {
                do_applyControl(discretisation, timeStep, rotorMeanWind, rotorTorque, rotorOmega, rotorWindVaneWind);
            }

            HOST_DEVICE_PREFIX void checkNeeds(bool & needsMeanVelocity, bool & needsTorque, bool & needsWindVene){
                do_checkNeeds(needsMeanVelocity, needsTorque, needsWindVene);
            }

            HOST_DEVICE_PREFIX walberla::mpi::MPISize getItemSize() {
                return do_getItemSize();
            }

            HOST_DEVICE_PREFIX void fillBuffer( BufferType_T * buffer ) {
                return do_fillBuffer(buffer);
            }

            HOST_DEVICE_PREFIX void unfillBuffer( BufferType_T * buffer ) {
                return do_unfillBuffer(buffer);
            }

        private:

            HOST_DEVICE_PREFIX virtual void do_print() const {
                printf("Base class ControlModel is not implemented.\n");
                exit(EXIT_FAILURE);
            }

            HOST_DEVICE_PREFIX virtual ControlModel * do_clone() const {
                printf("Base class ControlModel is not implemented.\n");
                exit(EXIT_FAILURE);
            }

            HOST_DEVICE_PREFIX virtual void do_applyControl(Discretisation * discretisation, const uint_t timeStep, 
                                                            const real_t rotorMeanWind, const real_t rotorTorque, 
                                                            const real_t rotorOmega, const Vector3<real_t> rotorWindVaneWind) {
                printf("Base class ControlModel is not implemented.\n");
                exit(EXIT_FAILURE);
            }

            HOST_DEVICE_PREFIX virtual void do_checkNeeds(bool & needsMeanVelocity, bool & needsTorque, bool & needsWindVane) {
                printf("Base class ControlModel is not implemented.\n");
                exit(EXIT_FAILURE);
            }

            HOST_DEVICE_PREFIX virtual walberla::mpi::MPISize do_getItemSize() {
                printf("Base class ControlModel is not implemented.\n");
                exit(EXIT_FAILURE);
            }

            HOST_DEVICE_PREFIX virtual void do_fillBuffer(BufferType_T * buffer) {
                printf("Base class ControlModel is not implemented.\n");
                exit(EXIT_FAILURE);
            }

            HOST_DEVICE_PREFIX virtual void do_unfillBuffer(BufferType_T * buffer) {
                printf("Base class ControlModel is not implemented.\n");
                exit(EXIT_FAILURE);
            }

        };
    }

    using control_model::ControlModel;

}

#endif //TURBINECORE_CONTROLMODEL_H
