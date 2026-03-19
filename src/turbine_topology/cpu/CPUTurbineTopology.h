
#pragma once

#ifndef TURBINECORE_CPUTURBINETOPOLOGY_H
#define TURBINECORE_CPUTURBINETOPOLOGY_H

#include "component/Component.h"
#include "domain/TurbineDomain.h"

#include "wind_turbine_core/ProjectDefines.h"

#include "walberla_helper/blockforest/BlockInfo.h"

#include "mpi/CustomMemoryBuffer.h"
#include "mpi/CPUMemoryAllocator.h"

#include <core/logging/Logging.h>
#include <core/mpi/MPIWrapper.h>
#include <core/mpi/Reduce.h>
#include <core/mpi/BufferSystem.h>
#include <domain_decomposition/IBlock.h>
#include <field/GhostLayerField.h>

#include "walberla_helper/field/interpolators/TrilinearFieldInterpolator.h"

namespace turbine_core {

    namespace topology {

        namespace cpu {

            template< typename T >
            class CPUTurbineTopology {

            public:

                using Buffer_T = mpi::CustomMemoryBuffer<mpi::CPUMemoryAllocator>;
                using BufferSystem_T = walberla::mpi::GenericBufferSystem<Buffer_T, Buffer_T>;
                using MPIInfo_T = std::map<walberla::mpi::MPIRank, walberla::mpi::MPISize>;

                HOST_DEVICE_PREFIX CPUTurbineTopology() {}
                HOST_DEVICE_PREFIX CPUTurbineTopology(const CPUTurbineTopology & ) = delete;
                HOST_DEVICE_PREFIX ~CPUTurbineTopology() {}

                HOST_PREFIX void callback( const Component::Function & function, const ComponentType & type, const uint_t timestep = 0 ) {
                    static_cast<T&>(*this).do_callback(function, type, timestep);
                }

                HOST_PREFIX void callback( const Component::Output & function, const walberla::filesystem::path & baseFolder,
                                           const uint_t timestep, const ComponentType & type ) {

                    WALBERLA_ROOT_SECTION() {

                        std::string filename;
                        switch (function) {
                            case Component::Output::GNUPLOT :
                                filename = "Gnuplot_" + std::to_string(timestep) + ".txt";
                                break;
                            case Component::Output::ORIENTATIONS :
                                filename = "Orientations_" + std::to_string(timestep) + ".txt";
                                break;
                            case Component::Output::FORCES : WALBERLA_ABORT("Not implemented.")
                        }

                        auto filepath = baseFolder / filename;
                        std::ofstream os(filepath.string(), std::ios::app);

                        if (!os.is_open()) WALBERLA_ABORT("Could not open file " << filepath.string() << ".")

                        static_cast<T &>(*this).do_callback(function, os);

                        os.close();
                    }
                }

                template< typename DensityField_T, typename VelocityField_T,
                          typename Interpolator_T>
                HOST_PREFIX void evaluateDensityAndVelocity( walberla::IBlock * block, std::shared_ptr<walberla::StructuredBlockForest> & storage,
                                                             const BlockDataID &densityFieldID,
                                                             const BlockDataID &velocityFieldID,
                                                             const uint_t &nPointsComponents,
                                                             const uint_t &nComponents) {


                    // get density and velocity fields
                    field::Field<real_t> densityField(block->getData<DensityField_T>(densityFieldID));
                    field::Field<real_t> velocityField(block->getData<VelocityField_T>(velocityFieldID));

                    /// EXTRACT NECESSARY DATA FROM BLOCK AND FOREST
                    blockforest::BlockInfo blockInfo{block, storage};

                    Interpolator_T densityInterpolator {blockInfo, &densityField};
                    Interpolator_T velocityInterpolator {blockInfo, &velocityField};

                    const uint_t nGL = densityField.nrOfGhostLayers();

                    int dummy{};
                    projectors::Projectors projectors(&densityInterpolator, &velocityInterpolator, &dummy, nGL);
                    for (uint_t pointIdy = 0; pointIdy < nPointsComponents; ++pointIdy) {
                        static_cast<T &>(*this).do_evaluateDensityAndVelocity(blockInfo, projectors, dummy, pointIdy);
                    }


                }

                HOST_PREFIX void applyControl( const uint_t & timeStep, const real_t & meanWindVelocity, const real_t & rotorTorque,
                    const real_t & rotorOmega, const Vector3<real_t> rotorWindVaneWind) {

                    static_cast<T&>(*this).do_applyControl(timeStep, meanWindVelocity, rotorTorque, rotorOmega, rotorWindVaneWind);

                }

                HOST_PREFIX void getControlNeeds(bool & needsMeanVelocity, bool & needsTorque, bool & needsWindVane){
                    static_cast<T&>(*this).do_getControlNeeds(needsMeanVelocity, needsTorque, needsWindVane);
                }

                HOST_PREFIX void syncNextNeighbour( const domain::TurbineDomain & domain, const mpi::SynchronisationMode mode,
                                                    const uint_t impactWidth, const uint_t turbineID, const bool ) {
                    // communication not needed for one MPI process
                    if(walberla::mpi::MPIManager::instance()->numProcesses() == 1)
                        return;

                    // current process does not contain turbine blocks
                    if( domain.getNumLocalAABBs() == 0 )
                        return;

                    // get receiver information for all processes
                    const auto & neighbourProcesses = domain.getNeighborProcesses();

                    if(neighbourProcesses.empty())
                        return;

                    BufferSystem_T bs( walberla::mpi::MPIManager::instance()->comm(), turbineID );

                    MPIInfo_T receiverInfo;
                    MPIInfo_T senderInfo;
                    std::map<walberla::mpi::MPIRank, std::vector<std::vector<uint_t>>> sendPoints;
                    std::map<walberla::mpi::MPIRank, std::vector<std::vector<uint_t>>> recvPoints;

                    for( uint_t nbProcessRank : neighbourProcesses ) {

                        const auto intRank = static_cast<int>(nbProcessRank);

                        auto & sendBuffer = bs.sendBuffer(nbProcessRank);
                        if (sendBuffer.isEmpty())
                        {
                            // fill empty buffers with a dummy byte to force transmission
                            receiverInfo[intRank] = walberla::mpi::MPISize(sizeof(walberla::uint8_t));
                            auto bufferPtr = sendBuffer.advance(1);
                            WALBERLA_ASSERT_NOT_NULLPTR(bufferPtr)
                            *bufferPtr = walberla::uint8_c(0);
                        }
                        static_cast<T&>(*this).do_setupCommunication(receiverInfo[intRank], senderInfo[intRank], nbProcessRank,
                                                                     sendPoints[intRank],
                                                                     recvPoints[intRank],
                                                                     domain, mode, impactWidth);

                        sendBuffer.resize(uint_t(senderInfo[intRank]) +
                                                              walberla::mpi::MPISize(sizeof(walberla::uint8_t)));

                    }

                    bs.setReceiverInfo(receiverInfo);

                    // fill buffer
                    for( uint_t nbProcessRank : neighbourProcesses ) {
                        auto & sendBuffer = bs.sendBuffer(nbProcessRank);
                        WALBERLA_ASSERT(!sendBuffer.isEmpty(), "Empty send buffer!" )
                        auto * pt = sendPoints[int(nbProcessRank)].data();
                        static_cast<T &>(*this).do_packTurbineData(sendBuffer, &pt, mode);
                    }

                    bs.sendAll();

                    // parse buffer data back to turbine
                    for( auto recv = bs.begin(); recv != bs.end(); ++recv ) {

                        walberla::uint8_t tmp;
                        tmp = *(recv.buffer().advanceNoResize(1));

                        auto * pt = recvPoints[int(recv.rank())].data();
                        static_cast<T &>(*this).do_unpackTurbineData(recv.buffer(), &pt, mode);
                    }
                }

                template< typename ForceField_T, typename ForceDistributor_T >
                HOST_PREFIX void spreadForces(walberla::IBlock * block, std::shared_ptr<walberla::StructuredBlockForest> & storage,
                                              const BlockDataID &forceFieldID,
                                              const uint_t &nPointsComponentsForSpread) {

                    /// EXTRACT NECESSARY DATA FROM BLOCK AND FOREST
                    blockforest::BlockInfo blockInfo{block, storage};

                    /// GET INTERPOLATORS AND DISTRIBUTOR
                    field::Field<real_t> forceField(block->getData<ForceField_T>(forceFieldID));
                    ForceDistributor_T forceDistributor {blockInfo, &forceField};

                    const uint_t nGL = forceField.nrOfGhostLayers();

                    int dummy{};
                    projectors::Projectors projectors(&dummy, &dummy, &forceDistributor, nGL);
                    for (uint_t pointIdy = 0; pointIdy < nPointsComponentsForSpread; pointIdy++) {
                        static_cast<T &>(*this).do_spreadForces(blockInfo, projectors, pointIdy);
                    }

                }

                template<typename Discretisation_T>
                HOST_PREFIX void addDiscretisation(const std::string & parentID, const std::string & componentID,
                                                   const ComponentType & type,
                                                   const std::shared_ptr<Discretisation_T> & discretisation)
                {
                    static_cast<T &>(*this).do_addDiscretisation(parentID, componentID, type, discretisation);
                }

                template<typename ForceModel_T, typename ControlModel_T>
                HOST_PREFIX void initialiseForceAndControlModel(const std::string & componentID,
                                                                const std::shared_ptr<ForceModel_T> & forceModel,
                                                                const std::shared_ptr<ControlModel_T> & controlModel)
                {
                    if(controlModel != nullptr) isControlled_ = true;
                    static_cast<T &>(*this).do_initialiseForceAndControlModel(componentID, forceModel, controlModel);
                }

                HOST_PREFIX void setupEnvironment(const std::string & baseName, const Vector3<real_t> & origin) {
                    static_cast<T&>(*this).do_setupEnvironment(baseName, origin);
                }

                HOST_PREFIX void deleteEnvironment() const {}

                HOST_PREFIX void writeForceOutput( walberla::IBlock * block, std::shared_ptr<walberla::StructuredBlockForest> &,
                                                   walberla::filesystem::path & baseFolder, const uint_t timestep ) {

                    static_cast<T*>(this)->do_writeForceOutput(baseFolder, timestep, block);

                }

                HOST_PREFIX void calculatePowerAndThrust(walberla::IBlock * block, std::shared_ptr<walberla::StructuredBlockForest> & storage,
                                                         real_t & power, real_t & thrust) {

                    blockforest::BlockInfo blockInfo{block, storage};

                    static_cast<T*>(this)->do_calculatePowerAndThrust(blockInfo, power, thrust);

                }

                HOST_PREFIX void calculateTorque(walberla::IBlock * block, std::shared_ptr<walberla::StructuredBlockForest> & storage,
                                                         real_t & torque) {

                    blockforest::BlockInfo blockInfo{block, storage};

                    static_cast<T*>(this)->do_calculateTorque(blockInfo, torque);

                }

                HOST_PREFIX void getOmega(real_t & omega) {
                    static_cast<T*>(this)->do_getHubOmega(omega);
                }

                HOST_PREFIX void calculateMeanVelocity(walberla::IBlock * block, std::shared_ptr<walberla::StructuredBlockForest> & storage,
                                                         real_t & totalVelocity, uint_t & nbPoints) {

                    blockforest::BlockInfo blockInfo{block, storage};

                    static_cast<T*>(this)->do_calculateMeanVelocity(blockInfo, totalVelocity, nbPoints);

                }

                HOST_PREFIX void calculateWindVaneVelocity(walberla::IBlock * block, std::shared_ptr<walberla::StructuredBlockForest> & storage,
                                         Vector3<real_t> & totalVelocity, uint_t & nbPoints) {

                    blockforest::BlockInfo blockInfo{block, storage};

                    static_cast<T*>(this)->do_calculateWindVaneVelocity(blockInfo, totalVelocity, nbPoints);

                }

                HOST_PREFIX AABB getAABB() {

                    Vector3<real_t> hubPosition{};
                    real_t bladeLength{};
                    real_t hubRadius{};
                    real_t towerHeight{};

                    static_cast<T*>(this)->do_getAABBInformation( towerHeight, hubPosition, hubRadius, bladeLength );

                    Vector3<real_t> halfWidth{bladeLength + hubRadius + towerHeight};

                    const auto min = hubPosition - halfWidth;
                    const auto max = hubPosition + halfWidth;

                    return {min, max};

                }

                HOST_PREFIX bool isControlled() const {return isControlled_;}
                HOST_PREFIX void setControlled(const bool isControlled) {isControlled_ = isControlled;}

            private:
                bool isControlled_{false};

            };

        }

    }

}

#endif //TURBINECORE_CPUTURBINETOPOLOGY_H
