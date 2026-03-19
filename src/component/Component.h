
#pragma once

#ifndef TURBINECORE_COMPONENT_H
#define TURBINECORE_COMPONENT_H

#include <memory>
#include <string>
#include <utility>

#include "ComponentType.h"
#include "discretisation/Discretisation.h"
#include "force_model/ForceModel.h"
#include "control_model/ControlModel.h"

#include "point3/Point3.h"

#include "wind_turbine_core/ProjectDefines.h"
#include "mpi/SynchronisationMode.h"

#include <core/logging/Logging.h>

namespace turbine_core {

    namespace blockforest {
        class BlockInfo;
    }

    namespace field {
        template<typename Type_T>
        class Field;
    }

    namespace projectors {
        class Projectors;
    }

    namespace component {

        HOST_DEVICE_PREFIX FORCEINLINE char * strClone(const char * str ) {
            size_t len = 0;
            do {} while (str[len++] != 0);
            ++len;
            auto tmp = new char [len]; // allocate for string and ending \0

            for(size_t i = 0; i < len; ++i) {
                tmp[i] = str[i];
            }

            return tmp;
        }

        class Component {

        public:

            enum Function {
                PRINT,
                UPDATE_DISCRETISATION,
                CALCULATE_FORCES,
                EVALUATE,
                RESET_COMMUNICATION_DATA
            };

            enum Output {
                GNUPLOT,
                ORIENTATIONS,
                FORCES
            };

            HOST_DEVICE_PREFIX Component( const char * ID, const ComponentType & type, Discretisation * discretisation,
                                          ForceModel * forceModel = nullptr, ControlModel * controlModel = nullptr,
                                          Point3<real_t> * reference = nullptr)
                    : componentID_(cloneID(ID)), componentType_(type),
                      referencePointFromParent_{reference}
            {
                if(discretisation) discretisation_ = discretisation->clone();
                if(forceModel) forceModel_ = forceModel->clone();
                if(controlModel) controlModel_ = controlModel->clone();
            }

            HOST_DEVICE_PREFIX Component( const Component & component )
                    : componentID_(cloneID(component.componentID_)), componentType_(component.componentType_),
                      discretisation_(component.discretisation_->clone()),
                      referencePointFromParent_(component.referencePointFromParent_) // non-owning -> do not deep copy!
            {
                if(component.forceModel_) forceModel_ = component.forceModel_->clone();
                if(component.controlModel_) controlModel_ = component.controlModel_->clone();
            }

            HOST_DEVICE_PREFIX ~Component() {
                delete[] componentID_;
                delete forceModel_;
                delete discretisation_;
                delete controlModel_;
            }

            HOST_DEVICE_PREFIX char * componentID() const {
                return &componentID_[0];
            }

            HOST_DEVICE_PREFIX ComponentType componentType() const {
                return componentType_;
            }

            HOST_DEVICE_PREFIX Discretisation * discretisation() const {
                return discretisation_;
            }

            HOST_DEVICE_PREFIX ForceModel * forceModel() const {
                return forceModel_;
            }

            HOST_DEVICE_PREFIX ControlModel * controlModel() const {
                return controlModel_;
            }

            HOST_DEVICE_PREFIX void initialiseForceAndControlModel( ForceModel * forceModel, ControlModel * controlModel ) {
                if(forceModel) forceModel_ = forceModel->clone();
                if(controlModel) controlModel_ = controlModel->clone();
            }

            HOST_DEVICE_PREFIX void setReferencePointFromParent( Point3<real_t> * referencePoint ) {
                referencePointFromParent_ = referencePoint;
            }

            HOST_DEVICE_PREFIX auto getReferencePointFromParent() {
                return referencePointFromParent_;
            }

            HOST_DEVICE_PREFIX auto getReferencePointForChild() {
                return discretisation_->getReferencePointForChild();
            }

            HOST_DEVICE_PREFIX void print() {
                printf("\nComponent %s (type = %i)\n", componentID_, componentType_);
                discretisation_->print();
                forceModel_->print();
                controlModel_->print();
            }

            HOST_DEVICE_PREFIX void updateDiscretisation(const uint_t timestep) {
                discretisation_->update(referencePointFromParent_, timestep);
            }

            HOST_DEVICE_PREFIX void applyControl(const uint_t timeStep, const real_t rotorMeanWind, const real_t rotorTorque,
                                                 const real_t rotorOmega, const Vector3<real_t> rotorWindVaneWind) {
                if(controlModel_) controlModel_->applyControl(discretisation_, timeStep, rotorMeanWind, rotorTorque, rotorOmega, rotorWindVaneWind);
            }

            HOST_DEVICE_PREFIX void getControlNeeds(bool & needsMeanVelocity, bool & needsTorque,
                bool & needsWindVane) {
                if(controlModel_) controlModel_->checkNeeds(needsMeanVelocity, needsTorque, needsWindVane);
            }

            HOST_DEVICE_PREFIX void updateForces() {
                if(forceModel_) forceModel_->calculateForces(discretisation_->points(), referencePointFromParent_);
            }

            HOST_DEVICE_PREFIX void spreadForces(const blockforest::BlockInfo & blockInfo,
                                                 const projectors::Projectors &projectors, uint_t pointIdy) {
                if (forceModel_) forceModel_->spreadForces(blockInfo, projectors, discretisation_->points(), pointIdy);
            }

            HOST_DEVICE_PREFIX void resetCommunicationData() {
                if(forceModel_) forceModel_->resetCommunicationData();
            }

            HOST_DEVICE_PREFIX void evaluateDensityAndVelocity( const blockforest::BlockInfo & blockInfo,
                                                                const projectors::Projectors &projectors,
                                                                uint_t pointIdx)
            {
                if (forceModel_)
                    forceModel_->evaluateDensityAndVelocity(blockInfo, projectors, discretisation_->points(), pointIdx);
            }

            HOST_DEVICE_PREFIX void callFunction( const Function & function, const ComponentType & type, const uint_t timestep = 0) {
                if(!(componentType_ & type))
                    return;

                switch (function) {

                    case Function::PRINT :
                        print();
                        break;
                    case Function::UPDATE_DISCRETISATION :
                        updateDiscretisation(timestep);
                        break;
                    case Function::CALCULATE_FORCES :
                        updateForces();
                        break;
                    case Function::RESET_COMMUNICATION_DATA :
                        resetCommunicationData();
                        break;
                    default:
                        printf("Invalid function for callback on Component.\n");
                        assert(0); // not implemented

                }
            }

            HOST_DEVICE_PREFIX void callFunction( const Output & function, const ComponentType & type ) {
                if(!(componentType_ & type))
                    return;

                switch (function) {

                    default:
                        printf("Invalid function for output callback on Component.\n");
                        assert(0); // not implemented

                }
            }

            template< typename Type_T >
            HOST_DEVICE_PREFIX void packData( Type_T * buffer, const mpi::SynchronisationMode mode, const uint_t idx ) {

                switch (mode) {
                    case mpi::SynchronisationMode::SYNC_ALL:
                        printf("SYNC_ALL is not yet implemented.");
                        break;
                    case mpi::SynchronisationMode::SYNC_FORCE: {
                        if(idx >= discretisation_->nPoints()) {
                            printf("idx %llu exceeds number of points.\n", idx);
                            assert(0);
                        }
                        if(forceModel_) {
                            buffer << forceModel_->points()[idx].forcesGlobal;
                        } else {
                            //FIXME waste of memory
                            buffer += sizeof(Vector3<real_t>);
                        }
                        break;
                    }
                    case mpi::SynchronisationMode::SYNC_MACRO: {
                        if(idx >= discretisation_->nPoints()) {
                            printf("idx %llu exceeds number of points.\n", idx);
                            assert(0);
                        }
                        if(forceModel_) {
                            auto &point = forceModel_->points()[idx];

                            memcpy(buffer, &point.density, sizeof(real_t));
                            buffer += sizeof(real_t);
                            buffer << point.velocityLocal;
                            buffer << point.windVelocityGlobal;
                        } else {
                            //FIXME waste of memory
                            buffer += 2 * sizeof(Vector3<real_t>);
                            buffer += sizeof(real_t);
                        }
                        break;
                    }
                    case mpi::SynchronisationMode::SYNC_GEOMETRY: {
                        if(idx >= discretisation_->nPoints()) {
                            printf("idx %llu exceeds number of points.\n", idx);
                            assert(0);
                        }
                        buffer << discretisation_->points()[idx];
                        break;
                    }
                    default:
                        printf("Invalid synchronisation mode for Component.\n");
                        assert(0); // not implemented
                }

            }

            template< typename Type_T >
            HOST_DEVICE_PREFIX void unpackData( Type_T * buffer, const mpi::SynchronisationMode mode, const uint_t idx ) {

                switch (mode) {
                    case mpi::SynchronisationMode::SYNC_ALL:
                        printf("SYNC_ALL is not yet implemented.");
                        break;
                    case mpi::SynchronisationMode::SYNC_FORCE: {
                        if(idx >= discretisation_->nPoints()) {
                            printf("idx %llu exceeds number of points.\n", idx);
                            assert(0);
                        }
                        if(forceModel_) {
                            buffer >> forceModel_->points()[idx].forcesGlobal;
                            if (forceModel_->remotePoint() != nullptr) {
                                forceModel_->remotePoint()[idx] = true;
                            }
                        } else {
                            //FIXME waste of memory
                            buffer += sizeof(Vector3<real_t>);
                        }
                        break;
                    }
                    case mpi::SynchronisationMode::SYNC_MACRO: {
                        if(idx >= discretisation_->nPoints()) {
                            printf("idx %llu exceeds number of points.\n", idx);
                            assert(0);
                        }

                        if(forceModel_) {
                            auto &point = forceModel_->points()[idx];

                            memcpy(&point.density, buffer, sizeof(real_t));
                            buffer += sizeof(real_t);

                            buffer >> point.velocityLocal;
                            buffer >> point.windVelocityGlobal;
                            if (forceModel_->remotePoint() != nullptr) {
                                forceModel_->remotePoint()[idx] = true;
                            }
                        } else {
                            //FIXME waste of memory
                            buffer += 2 * sizeof(Vector3<real_t>);
                            buffer += sizeof(real_t);
                        }
                        break;
                    }
                    case mpi::SynchronisationMode::SYNC_GEOMETRY: {
                        if(idx >= discretisation_->nPoints()) {
                            printf("idx %llu exceeds number of points.\n", idx);
                            assert(0);
                        }
                        buffer >> discretisation_->points()[idx];
                        break;
                    }
                    default:
                        printf("Invalid synchronisation mode for Component.\n");
                        assert(0); // not implemented
                }

            }

            template< typename Type_T >
            HOST_DEVICE_PREFIX void packControlData( Type_T * buffer ) {
                if(controlModel_) {
                    controlModel_->fillBuffer(buffer);
                }
            }


            template< typename Type_T >
            HOST_DEVICE_PREFIX void unpackControlData( Type_T * buffer ) {
                if(controlModel_) {
                    controlModel_->unfillBuffer(buffer);
                }
            }

            HOST_DEVICE_PREFIX bool isControlled() {
                return controlModel_ != nullptr;
            }

            HOST_DEVICE_PREFIX walberla::mpi::MPISize getControlSize() {
                return controlModel_->getItemSize();
            }

        private:

            HOST_DEVICE_PREFIX static char * cloneID(const char * str ) {
                size_t len = 0;
                do {} while (str[len++] != 0);
                ++len; // account for terminating char \0
                auto tmp = new char [len];

                // do not copy null terminator...
                for(size_t i = 0; i < len - 1; ++i) {
                    tmp[i] = str[i];
                }

                return tmp;
            }

            /// MEMBER VARIABLES

            char * componentID_{nullptr};

            ComponentType componentType_;

            ForceModel * forceModel_{nullptr};
            Discretisation * discretisation_{nullptr};
            ControlModel * controlModel_{nullptr};

            // non-owing!
            Point3<real_t> * referencePointFromParent_{};

        };

    }

    using component::Component;

}

#endif //TURBINECORE_COMPONENT_H
