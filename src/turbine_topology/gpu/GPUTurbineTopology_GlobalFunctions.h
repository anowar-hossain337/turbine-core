#pragma once

#ifndef TURBINECORE_GPUTURBINETOPOLOGY_GLOBALFUNCTIONS_H
#define TURBINECORE_GPUTURBINETOPOLOGY_GLOBALFUNCTIONS_H

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

#include "component/ComponentType.h"
#include "component/Component.h"

#include "walberla_helper/field/Field.h"
#include "walberla_helper/field/Projectors.h"

#include "wind_turbine_core/math/Vector3.h"
#include "discretisation/Disk.h"

namespace turbine_core {

    namespace topology {

        namespace gpu {

            namespace internal {
                MANAGED_PREFIX DEVICE_PREFIX uint_t nComponents{0};
                MANAGED_PREFIX DEVICE_PREFIX uint_t nControlledComponents{0};
                DEVICE_PREFIX uint_t * nPoints{nullptr};
                MANAGED_PREFIX DEVICE_PREFIX uint_t nBlades{0};
                MANAGED_PREFIX DEVICE_PREFIX uint_t nBladePoints{0};
            }

            template<typename T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void setupOnGPU(T ** base, const char * baseName, const Vector3<real_t> origin) {
                auto point = new discretisation::Disk(uint_t(0), origin, Vector3<real_t>(), Quaternion<real_t>());
                ForceModel * noneForceModel = nullptr;
                ControlModel * noneControlModel = nullptr;
                (*base) = new T(baseName, ComponentType::BASE, nullptr, point, noneForceModel, noneControlModel);
                delete point;
                delete noneForceModel;
                delete noneControlModel;
            }

            template< typename T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void deleteOnGPU(T ** base) {
                delete (*base);
            }

            template< typename T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void globalCallbackFct( T ** RESTRICT base, const Component::Function function, const ComponentType type, const uint_t timestep ) {
                (*base)->do_callback(function, type, timestep);
            }

            template< typename T, typename DataType_T >
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void getOutputDataOnGPU( T ** RESTRICT base, const Component::Output function, DataType_T * data ) {
                (*base)->do_getOutputData(function, data);
            }

            template< typename T >
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void getControlNeedsOnGPU( T ** RESTRICT base, bool * needsMeanVelocity, bool * needsTorque,
            bool * needsWindVane) {
                (*base)->do_getControlNeeds( *needsMeanVelocity, *needsTorque, *needsWindVane );
            }

            template< typename T >
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void applyControlOnGPU( T ** RESTRICT base, const uint_t timeStep,
                                                  real_t rotorMeanVelocity, real_t rotorTorque, real_t rotorOmega,
												  Vector3<real_t> rotorWindVaneWind) {
                (*base)->do_applyControl( timeStep, rotorMeanVelocity, rotorTorque, rotorOmega, rotorWindVaneWind );
            }

            template< typename T, typename Interpolator_T >

            GLOBAL_PREFIX void evaluateDensityAndVelocityOnGPU( T ** RESTRICT base, const blockforest::BlockInfo blockInfo,
                                                                field::Field <real_t> densityField,
                                                                field::Field <real_t> velocityField,
                                                                const uint_t nComponentsBlade) {
                Interpolator_T densityInterpolator {blockInfo, &densityField};
                Interpolator_T velocityInterpolator {blockInfo, &velocityField};

                const uint_t nGL = densityField.nrOfGhostLayers();

                int dummy;
                projectors::Projectors projectors{ &densityInterpolator, &velocityInterpolator, &dummy, nGL };

                int pointIdx = blockIdx.x;
                int pointIdy = threadIdx.x;

                if (pointIdy < nComponentsBlade) {
                    (*base)->do_evaluateDensityAndVelocity(blockInfo, projectors, pointIdx, pointIdy);
                }


            }

            template<typename T, typename ForceDistributor_T>

            GLOBAL_PREFIX void spreadForcesOnGPU( T ** RESTRICT base, const blockforest::BlockInfo blockInfo,
                                                  field::Field <real_t> forceField,
                                                  const uint_t nPointsComponentsForSpread) {
                ForceDistributor_T forceDistributor{ blockInfo, &forceField };

                const uint_t nGL = forceField.nrOfGhostLayers();

                int dummy{};
                uint_t pointIdy = blockIdx.x * blockDim.x + threadIdx.x;

                projectors::Projectors projectors(&dummy, &dummy, &forceDistributor, nGL);
                if (pointIdy < nPointsComponentsForSpread)
                    (*base)->do_spreadForces(blockInfo, projectors, pointIdy);
            }

            template< typename Discretisation_T, typename T >
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void addGPUDiscretisation(T ** RESTRICT base,
                                                    const char * RESTRICT parentID, const char * RESTRICT componentID,
                                                    const ComponentType type,
                                                    Discretisation_T const*const* RESTRICT p_discretisation ) {

                //NOTE necessary because as virtual functions cannot be handled for p_* (keeps stuck in infinite loop...)
                auto discretisation = new Discretisation_T(**p_discretisation);

                (*base)->do_addDiscretisation(parentID, componentID, type, discretisation);

                delete discretisation;

            }

            template< typename ForceModel_T, typename ControlModel_T, typename T >
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void initialiseGPUForceAndControlModel(T ** RESTRICT base,
                                                                 const char * RESTRICT componentID,
                                                                 const bool provideForceModel,
                                                                 ForceModel_T const*const* RESTRICT p_forceModel,
                                                                 const bool provideControlModel,
                                                                 ControlModel_T const*const* RESTRICT p_controlModel ) {

                ForceModel_T * forceModel{nullptr};
                ControlModel_T * controlModel{nullptr};

                //NOTE must provide bools as variables as check via nullptr is not possible
                if(provideForceModel) forceModel = new ForceModel_T(**p_forceModel);
                if(provideControlModel) controlModel = new ControlModel_T(**p_controlModel);

                (*base)->do_initialiseForceAndControlModel(componentID, forceModel, controlModel);

                delete forceModel;
                delete controlModel;
            }

            template< typename T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void setupCommunicationOnGPU(T ** base,
                                                       const uint_t ownRank, const uint_t neighborRank,
                                                       int * sendPoints, int * recvPoints,
                                                       const uint_t nLocalAABBs, const uint_t nNeighborAABBs,
                                                       math::AABB * const localAABBs, math::AABB * const neighborAABBs,
                                                       uint_t * nSendPoints, walberla::mpi::MPISize * receiverSize,
                                                       const mpi::SynchronisationMode mode, const real_t impactWidth) {
                (*base)->do_setupCommunication(
                        ownRank, neighborRank,
                        sendPoints, recvPoints,
                        nLocalAABBs, nNeighborAABBs,
                        localAABBs, neighborAABBs,
                        *nSendPoints, *receiverSize,
                        mode, impactWidth);
            }

            template< typename T, typename Type_T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void packTurbineDataOnGPU(T ** base,
                                                    Type_T * buffer,
                                                    int const * const pt,
                                                    const mpi::SynchronisationMode mode) {

                (*base)->do_packTurbineData(buffer, pt, mode);

            }

            template< typename T, typename Type_T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void unpackTurbineDataOnGPU(T ** base,
                                                      Type_T * buffer,
                                                      int * pt,
                                                      const mpi::SynchronisationMode mode) {

                (*base)->do_unpackTurbineData(buffer, pt, mode);

            }

            template< typename T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void getTurbineCharacteristicsOnGPU(T ** base, uint_t * pointPtr) {
                internal::nComponents = (*base)->do_getNComponents();
                internal::nControlledComponents = (*base)->do_getNControlledComponents();

                (*base)->do_getNPoints( pointPtr );
            }

            template< typename T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void getBladeDataOnGPU(T ** base) {
                (*base)->do_getBladeData();
            }

            template< typename T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void getForceDataOnGPU(T ** base, const blockforest::BlockInfo blockInfo,
                                                 uint_t * const nLocalPoints, uint_t * const localPoints,
                                                 real_t * const azimuth, Vector2<real_t> * const localForces,
                                                 Vector3<real_t> * const localAero, real_t * const relVelocity) {
                (*base)->do_getForceData(blockInfo, nLocalPoints, localPoints, azimuth, localForces, localAero, relVelocity);
            }

            template< typename T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void calculatePowerAndThrustOnGPU(T ** base, const blockforest::BlockInfo blockInfo, real_t * power, real_t * thrust) {
                (*base)->do_calculatePowerAndThrust(blockInfo, *power, *thrust);
            }

            template< typename T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void calculateTorqueOnGPU(T ** base, const blockforest::BlockInfo blockInfo, real_t * torque) {
                (*base)->do_calculateTorque(blockInfo, *torque);
            }

            template< typename T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void getOmegaOnGPU(T ** base, real_t * omega) {
                (*base)->do_getOmega(*omega);
            }

            template< typename T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void calculateMeanVelocityOnGPU(T ** base, const blockforest::BlockInfo blockInfo, real_t * meanVelocity, uint_t * nPoints) {
                (*base)->do_calculateMeanVelocity(blockInfo, *meanVelocity, *nPoints);
            }

            template< typename T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void calculateWindVaneVelocityOnGPU(T ** base, const blockforest::BlockInfo blockInfo, Vector3<real_t> * windVaneVelocity, uint_t * nPoints) {
                (*base)->do_calculateWindVaneVelocity(blockInfo, *windVaneVelocity, *nPoints);
            }

            template< typename T>
            LAUNCH_BOUNDS(1)
            GLOBAL_PREFIX void getAABBInformationOnGPU(T ** base, real_t * info) {
                (*base)->do_getAABBInformation(info);
            }

        } // namespace gpu

    } // namespace data_structure

} // namespace turbine_core

#endif //TURBINECORE_GPUTURBINETOPOLOGY_GLOBALFUNCTIONS_H
