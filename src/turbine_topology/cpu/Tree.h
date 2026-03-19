#pragma once

#ifndef TURBINECORE_CPU_TREE_H
#define TURBINECORE_CPU_TREE_H

#include <core/Filesystem.h>

#include "turbine_topology/cpu/CPUTurbineTopology.h"

#include "discretisation/Disk.h"
#include "discretisation/Line.h"
#include "control_model/RAWSOmegaModel.h"

#include <vector>
#include "conversion/Conversion.h"
#include "force_model/ActuatorLineModel.h"
#include "mpi/SynchronisationMode.h"

#include "domain/TurbineDomain.h"


namespace turbine_core {

    namespace topology {

        namespace cpu {

            class Tree final : public CPUTurbineTopology<Tree> {

                friend class CPUTurbineTopology<Tree>;

            public:

                HOST_PREFIX Tree() {}

                HOST_PREFIX Tree( const char * ID, const ComponentType & type, Tree * parent,
                                  discretisation::Discretisation * discretisation )
                        : component{std::make_shared<component::Component>(ID, type, discretisation)}, parent_(parent)
                {
                    if(parent_ != nullptr) {
                        component->setReferencePointFromParent(parent_->component->getReferencePointForChild());
                    }
                }

                HOST_PREFIX Tree( const Tree & tree )
                        : component(tree.component), parent_(tree.parent_), children_(tree.children_)
                {}

            private:

                HOST_PREFIX void do_setupEnvironment(const std::string & baseName, const Vector3<real_t> & origin) {
                    auto point = std::make_shared<discretisation::Disk>(real_t(0), origin);
                    std::shared_ptr<force_model::ForceModel> forceModel = nullptr;
                    std::shared_ptr<control_model::ControlModel> controlModel = nullptr;

                    component = std::make_shared<Component>(baseName.c_str(), ComponentType::BASE, point.get(), forceModel.get(), controlModel.get());
                }

                HOST_PREFIX void do_callback( const Component::Function & function, const ComponentType & type) const {

                    if(!(component->componentType() & ComponentType::BASE))
                        component->callFunction(function, type);

                    for(auto& child : children_) {
                        child.do_callback(function, type);
                    }
                }

                HOST_PREFIX void do_callback( const Component::Function & function, const ComponentType & type, const uint_t timestep) const {

                    if(!(component->componentType() & ComponentType::BASE))
                        component->callFunction(function, type, timestep);

                    for(auto& child : children_) {
                        child.do_callback(function, type, timestep);
                    }
                }

                HOST_PREFIX void do_callback( const Component::Output & function, std::ofstream & os ) const {

                    auto discretisation = component->discretisation();
                    auto points = discretisation->points();

                    for(uint_t p = 0; p < discretisation->nPoints(); ++p) {

                        if( (p % 3 != 0) && p != discretisation->nPoints() - 1 )
                            continue;

                        auto position = points[p].position;
                        auto matrix = points[p].orientation.toRotationMatrix();
                        auto dx = Vector3<real_t>(matrix[0], matrix[3], matrix[6]);
                        auto dy = Vector3<real_t>(matrix[1], matrix[4], matrix[7]);
                        auto dz = Vector3<real_t>(matrix[2], matrix[5], matrix[8]);

                        switch (function) {
                            case Component::Output::GNUPLOT :
                                os << position[0] << "\t" << position[1] << "\t" << position[2] << "\n";
                                break;
                            case Component::Output::ORIENTATIONS :

                                os << position[0] << "\t" << position[1] << "\t" << position[2] << "\t"
                                   << dx[0] << "\t" << dx[1] << "\t" << dx[2] << "\tx\n";
                                os << position[0] << "\t" << position[1] << "\t" << position[2] << "\t"
                                   << dy[0] << "\t" << dy[1] << "\t" << dy[2] << "\ty\n";
                                os << position[0] << "\t" << position[1] << "\t" << position[2] << "\t"
                                   << dz[0] << "\t" << dz[1] << "\t" << dz[2] << "\tz\n";
                                break;
                            default :
                                break;
                        }
                    }

                    os << "\n\n";

                    for(auto& child : children_) {
                        child.do_callback(function, os);
                    }
                }

                HOST_PREFIX void do_addDiscretisation (const std::string & parentID, const std::string & ID, const ComponentType & type,
                                                       const std::shared_ptr<discretisation::Discretisation> & discretisation ) {

                    const std::string componentIDstr{component->componentID()};

                    if(parentID == componentIDstr) {
                        children_.emplace_back(ID.c_str(), type, this, discretisation.get());
                    } else {
                        for(auto& child : children_) {
                            child.do_addDiscretisation(parentID, ID, type, discretisation);
                        }
                    }

                    //TODO how to check if child could be added?

                }

                HOST_PREFIX void do_initialiseForceAndControlModel (const std::string & ID,
                                                                    const std::shared_ptr<force_model::ForceModel> & forceModel,
                                                                    const std::shared_ptr<control_model::ControlModel> & controlModel) {

                    const std::string componentIDstr{component->componentID()};

                    ForceModel * fm{nullptr};
                    ControlModel * cm{nullptr};

                    if(forceModel) fm = forceModel.get();
                    if(controlModel) cm = controlModel.get();

                    if(ID == componentIDstr) {
                        component->initialiseForceAndControlModel(fm, cm);
                    } else {
                        for(auto& child : children_) {
                            child.do_initialiseForceAndControlModel(ID, forceModel, controlModel);
                        }
                    }

                }

                HOST_PREFIX void do_print() const {

                    printf("\nComponent: %s (type = %iu)\n", component->componentID(), component->componentType());
                    component->print();

                    for(auto& child : children_) {
                        child.do_print();
                    }
                }

                HOST_PREFIX void do_applyControl(const uint_t & timeStep, const real_t rotorMeanWind, const real_t rotorTorque, const real_t rotorOmega, const Vector3<real_t> rotorWindVaneWind)
                {
                    if(component->isControlled()) {
                        component->applyControl(timeStep, rotorMeanWind, rotorTorque, rotorOmega, rotorWindVaneWind);
                    }

                    for (auto &child: children_) {
                        child.do_applyControl(timeStep, rotorMeanWind, rotorTorque, rotorOmega, rotorWindVaneWind);
                    }
                }

                HOST_PREFIX void do_getControlNeeds(bool & needsMeanVelocity, bool & needsTorque, bool & needsWindVane)
                {
                    if(component->isControlled()) {
                        bool meanVel{false}, torque{false}, windVane{false};
                        component->getControlNeeds(meanVel, torque, windVane);
                        needsMeanVelocity += meanVel;
                        needsTorque += torque;
                        needsWindVane += windVane;
                    }

                    for (auto &child: children_) {
                        child.do_getControlNeeds(needsMeanVelocity, needsTorque, needsWindVane);
                    }
                }

                HOST_PREFIX void do_evaluateDensityAndVelocity( const blockforest::BlockInfo & blockInfo,
                                                                const projectors::Projectors &projectors,
                                                                uint_t pointIdx,
                                                                uint_t pointIdy)
                {
                    component->evaluateDensityAndVelocity(blockInfo, projectors, pointIdy);

                    for(auto& child : children_) {
                        child.do_evaluateDensityAndVelocity(blockInfo, projectors, pointIdx, pointIdy);
                    }
                }

                HOST_PREFIX void do_spreadForces( const blockforest::BlockInfo & blockInfo,
                                                  const projectors::Projectors &projectors,
                                                  uint_t pointIdy)
                {
                    component->spreadForces(blockInfo, projectors, pointIdy);

                    for(auto& child : children_) {
                        child.do_spreadForces(blockInfo, projectors, pointIdy);
                    }
                }

                HOST_PREFIX void do_setupCommunication(walberla::mpi::MPISize & receiverMPISize, walberla::mpi::MPISize & senderMPISize,
                                                       const uint_t process,
                                                       std::vector<std::vector<uint_t>> & sendPoints,
                                                       std::vector<std::vector<uint_t>> & recvPoints,
                                                       const domain::TurbineDomain & domain,
                                                       const mpi::SynchronisationMode mode, const uint_t impactWidth) {

                    sendPoints.emplace_back();
                    recvPoints.emplace_back();

                    const auto ownRank = static_cast<uint_t>(walberla::mpi::MPIManager::instance()->rank());

                    auto * discretisation = component->discretisation();
                    WALBERLA_ASSERT_NOT_NULLPTR(discretisation, "Tried to access a blade with invalid discretisation!")

                    const uint_t nPoints = discretisation->nPoints();
                    const auto * points = discretisation->points();

                    // loop over all points to get data size to communicate
                    for( uint_t p = 0; p < nPoints; ++p ) {
                        const auto pos = walberla::Vector3<real_t>(points[p].position.data());

                        // point lies in local block but impact width reaches into remote block
                        //   -> send to remote neighbour
                        if( domain.intersectsWithProcessSubdomain(ownRank, pos, real_t(0)) &&
                            domain.intersectsWithProcessSubdomain(process, pos, real_t(impactWidth)) ) {
                            sendPoints.back().emplace_back(p);

                            senderMPISize += walberla::mpi::MPISize(mpi::itemSize(mode));
                            if(p == 0 && component->isControlled()){
                                //FIXME currently done in every synchronisation, not only for SYNC_CONTROL
                                senderMPISize += walberla::mpi::MPISize(component->getControlSize());
                            }
                        }

                        // point lies in remote block but impact width reaches into local block
                        //   -> receive from neighbour
                        if( domain.intersectsWithProcessSubdomain(ownRank, pos, real_t(impactWidth)) &&
                            domain.intersectsWithProcessSubdomain(process, pos, real_t(0)) ) {
                            recvPoints.back().emplace_back(p);

                            receiverMPISize += walberla::mpi::MPISize(mpi::itemSize(mode));
                            if(p == 0 && component->isControlled()){
                                //FIXME currently done in every synchronisation, not only for SYNC_CONTROL
                                receiverMPISize += walberla::mpi::MPISize(component->getControlSize());
                            }
                        }
                    }

                    for( auto & child : children_ ) {
                        child.do_setupCommunication(receiverMPISize, senderMPISize, process, sendPoints, recvPoints, domain, mode, impactWidth);
                    }

                }

                HOST_PREFIX void do_packTurbineData(Buffer_T & buffer, std::vector<uint_t> ** sendPoints,
                                                    const mpi::SynchronisationMode mode) {

                    const auto size = static_cast<size_t>(mpi::itemSize(mode));

                    for( uint_t idx = 0; idx < (*sendPoints)->size(); ++idx ) {
                        WALBERLA_ASSERT_GREATER_EQUAL(buffer.allocSize() - buffer.size(), size, "Not sufficient space in buffer")
                        auto bufferPtr = buffer.advanceNoResize(size);
                        component->packData(bufferPtr, mode, (*sendPoints)->at(idx));

                        if((*sendPoints)->at(idx) == 0 && component->isControlled()){
                            //FIXME currently done in every synchronisation, not only for SYNC_CONTROL
                            const auto controlSize = static_cast<size_t>(component->getControlSize());
                            WALBERLA_ASSERT_GREATER_EQUAL(buffer.allocSize() - buffer.size(), controlSize, "Not sufficient space in buffer")
                            auto bufferPtr = buffer.advanceNoResize(controlSize);
                            component->packControlData(bufferPtr);
                        }
                    }
                    ++(*sendPoints);

                    for( auto & child : children_ ) {
                        child.do_packTurbineData(buffer, sendPoints, mode);
                    }

                }

                HOST_PREFIX void do_unpackTurbineData(Buffer_T & buffer, std::vector<uint_t> ** recvPoints,
                                                      const mpi::SynchronisationMode mode) {

                        const auto size = static_cast<size_t>(mpi::itemSize(mode));

                        for( uint_t idx = 0; idx < (*recvPoints)->size(); ++idx ) {

                            WALBERLA_ASSERT_GREATER_EQUAL(buffer.allocSize() - buffer.size(), size, "Not sufficient data in buffer")
                            auto bufferPtr = buffer.advanceNoResize(size);
                            component->unpackData(bufferPtr, mode, (*recvPoints)->at(idx));

                            if((*recvPoints)->at(idx) == 0 && component->isControlled()){
                                //FIXME currently done in every synchronisation, not only for SYNC_CONTROL
                                const auto controlSize = static_cast<size_t>(component->getControlSize());
                                WALBERLA_ASSERT_GREATER_EQUAL(buffer.allocSize() - buffer.size(), controlSize, "Not sufficient space in buffer")
                                auto bufferPtr = buffer.advanceNoResize(controlSize);
                                component->unpackControlData(bufferPtr);
                            }
                        }

                        ++(*recvPoints);

                    for( auto & child : children_ ) {
                        child.do_unpackTurbineData(buffer, recvPoints, mode);
                    }

                }

                HOST_PREFIX void do_writeForceOutput( const walberla::filesystem::path & baseFolder, const uint_t timestep, walberla::IBlock * block ) {

                    if(!!(component->componentType() & ComponentType::BLADE)) {

                        auto bladeFolder = baseFolder / walberla::filesystem::path(component->componentID());

                        auto * line = dynamic_cast<discretisation::Line*>(component->discretisation());
                        WALBERLA_ASSERT_NOT_NULLPTR(line, "Tried to access a blade with non-Line discretisation!")

                        auto * reference = parent_->component->getReferencePointFromParent();

                        auto * discretisationPoints = component->discretisation()->points();
                        auto * forcePoints = component->forceModel()->points();

                        const real_t Cf = Conversion::C_m() / ( Conversion::C_t() * Conversion::C_t() );
                        const real_t factor = Cf / line->elementLength_;

                        for( uint_t p = 0; p < line->nPoints(); ++p ) {

                            auto & position = discretisationPoints[p].position;

                            if(!block->getAABB().contains(position[0], position[1], position[2]))
                                continue;

                            const std::string pString{std::to_string(p)};

                            const std::string filename{ "Point_" + std::string(3 - pString.length(), '0') + pString + ".txt" };

                            auto filepath = bladeFolder / filename;

                            std::ofstream os ( filepath.string(), std::ios::app );

                            if(os.is_open()) {

                                real_t azimuthAngle{};
                                line->calculateAzimuthAngle(reference, p, azimuthAngle);

                                os << timestep << "\t";
                                os << azimuthAngle << "\t";

                                Vector3<real_t> localForces = discretisationPoints[p].orientation.getInverse().rotate(forcePoints[p].forcesGlobal);

                                auto aoa = forcePoints[p].angleOfAttack;
                                auto polars = forcePoints[p].polar(aoa);

                                auto relativeVelocity = forcePoints[p].velocityLocal;
                                relativeVelocity[2] = 0.;
                                real_t velocityMagnitude = relativeVelocity.length();

                                os << localForces[0] * factor << "\t" << localForces[1] * factor << "\t" << forcePoints[p].angleOfAttack / walberla::math::pi * real_t(180)
                                   << "\t" << polars[0] << "\t" << polars[1] << "\t" << velocityMagnitude * Conversion::C_l() / Conversion::C_t();

                                os << "\n";

                                os.close();
                            } else {
                                WALBERLA_ABORT("Could not open file " + filepath.string())
                            }

                        }
                    }

                    for( auto & child : children_ ) {
                        child.do_writeForceOutput(baseFolder, timestep, block);
                    }

                }

                HOST_PREFIX void do_calculateMeanVelocity( const blockforest::BlockInfo & blockInfo, real_t & totalVelocity, uint_t & nbPoints ){

                    if(!!(component->componentType() & ComponentType::BLADE)) {

                        WALBERLA_ASSERT_NOT_NULLPTR(dynamic_cast<discretisation::Line*>(component->discretisation()), "Tried to access a blade with non-Line discretisation!")

                        auto * discretisationPoints = component->discretisation()->points();
                        WALBERLA_ASSERT_NOT_NULLPTR(discretisationPoints, "Component returned null pointer as discretisation points.")

                        auto * forcePoints = component->forceModel()->points();
                        WALBERLA_ASSERT_NOT_NULLPTR(forcePoints, "Component returned null pointer as actuator data points.")

                        const uint_t nPoints = component->discretisation()->nPoints();

                        auto reference = component->getReferencePointFromParent();

                        auto orientation = reference->orientation;

                        const auto & blockBB = blockInfo.getAABB();

                        for(uint_t p = 0; p < nPoints; ++p) {
                            const auto position = discretisationPoints[p].position;
                            if(blockBB.contains(position)) {
                                ++nbPoints;
                                // totalVelocity is local
                                totalVelocity -= orientation.getInverse().rotate(forcePoints[p].windVelocityGlobal)[2]; // Velocity along the hub axis
                            }
                        }

                    } // if blade

                    for( auto & child : children_ ) {
                        child.do_calculateMeanVelocity(blockInfo, totalVelocity, nbPoints);
                    }


                } // function do_calculateMeanVelocity

                HOST_PREFIX void do_calculateWindVaneVelocity( const blockforest::BlockInfo & blockInfo, Vector3<real_t> & totalVelocity, uint_t & nbPoints ){

                    if(!!(component->componentType() & ComponentType::BLADE)) {

                        WALBERLA_ASSERT_NOT_NULLPTR(dynamic_cast<discretisation::Line*>(component->discretisation()), "Tried to access a blade with non-Line discretisation!")

                        auto * discretisationPoints = component->discretisation()->points();
                        WALBERLA_ASSERT_NOT_NULLPTR(discretisationPoints, "Component returned null pointer as discretisation points.")

                        auto * forcePoints = component->forceModel()->points();
                        WALBERLA_ASSERT_NOT_NULLPTR(forcePoints, "Component returned null pointer as actuator data points.")

                        const uint_t nPoints = component->discretisation()->nPoints();

                        auto reference = component->getReferencePointFromParent();

                        auto orientation = reference->orientation;

                        const auto & blockBB = blockInfo.getAABB();

                        const auto position = discretisationPoints[0].position;
                        if(blockBB.contains(position)){
                            ++nbPoints;
                            // totalVelocity is global
                            totalVelocity += forcePoints[0].windVelocityGlobal; // Velocity in the global reference frame
                        }
                    } // if blade

                    for( auto & child : children_ ) {
                        child.do_calculateWindVaneVelocity(blockInfo, totalVelocity, nbPoints);
                    }
                } // function do_calculateWindVaneVelocity

                HOST_PREFIX void do_calculatePowerAndThrust( const blockforest::BlockInfo & blockInfo, real_t & power, real_t & thrust ) {

                    if(!!(component->componentType() & ComponentType::BLADE)) {

                        WALBERLA_ASSERT_NOT_NULLPTR(dynamic_cast<discretisation::Line*>(component->discretisation()), "Tried to access a blade with non-Line discretisation!")

                        auto * discretisationPoints = component->discretisation()->points();
                        WALBERLA_ASSERT_NOT_NULLPTR(discretisationPoints, "Component returned null pointer as discretisation points.")

                        auto * forcePoints = component->forceModel()->points();
                        WALBERLA_ASSERT_NOT_NULLPTR(forcePoints, "Component returned null pointer as actuator data points.")

                        const uint_t nPoints = component->discretisation()->nPoints();
                        auto reference = component->getReferencePointFromParent();

                        auto orientation = reference->orientation;

                        auto hub = dynamic_cast<discretisation::Disk*>(this->parent_->component->discretisation());
                        auto rotVelocity = hub->relativeRotationalVelocity();

                        const real_t Cf = Conversion::C_m() * Conversion::C_l() / ( Conversion::C_t() * Conversion::C_t() );
                        const real_t thrustFactor = -Cf;
                        const real_t powerFactor = +Cf * Conversion::C_l() / Conversion::C_t() * rotVelocity[2];

                        const auto & blockBB = blockInfo.getAABB();

                        for(uint_t p = 0; p < nPoints; ++p) {
                            const auto position = discretisationPoints[p].position;
                            if(blockBB.contains(position)) {
                                // thrust and power are local
                                thrust += orientation.getInverse().rotate(forcePoints[p].forcesGlobal)[2] * thrustFactor;
                                power += orientation.getInverse().rotate((position - reference->position) % forcePoints[p].forcesGlobal)[2] * powerFactor;
                            }
                        }

                    } // if blade

                    for( auto & child : children_ ) {
                        child.do_calculatePowerAndThrust(blockInfo, power, thrust);
                    }
                } // function do_calculatePowerAndThrust

                HOST_PREFIX void do_calculateTorque( const blockforest::BlockInfo & blockInfo, real_t & torque ) {

                    if(!!(component->componentType() & ComponentType::BLADE)) {

                        WALBERLA_ASSERT_NOT_NULLPTR(dynamic_cast<discretisation::Line*>(component->discretisation()), "Tried to access a blade with non-Line discretisation!")

                        auto * discretisationPoints = component->discretisation()->points();
                        WALBERLA_ASSERT_NOT_NULLPTR(discretisationPoints, "Component returned null pointer as discretisation points.")

                        auto * forcePoints = component->forceModel()->points();
                        WALBERLA_ASSERT_NOT_NULLPTR(forcePoints, "Component returned null pointer as actuator data points.")

                        const uint_t nPoints = component->discretisation()->nPoints();
                        auto reference = component->getReferencePointFromParent();

                        auto orientation = reference->orientation;

                        auto hub = dynamic_cast<discretisation::Disk*>(this->parent_->component->discretisation());
                        auto rotVelocity = hub->relativeRotationalVelocity();

                        const auto & blockBB = blockInfo.getAABB();

                        for(uint_t p = 0; p < nPoints; ++p) {
                            const auto position = discretisationPoints[p].position;
                            if(blockBB.contains(position)) {
                                // torque is local
                                torque += orientation.getInverse().rotate((position - reference->position) % forcePoints[p].forcesGlobal)[2];
                            }
                        }

                    } // if blade

                    for( auto & child : children_ ) {
                        child.do_calculateTorque(blockInfo, torque);
                    }
                } // function do_calculateTorque

                HOST_PREFIX void do_getHubOmega( real_t & omega ) {

                    if(!!(component->componentType() & ComponentType::HUB)) {

                        WALBERLA_ASSERT_NOT_NULLPTR(dynamic_cast<discretisation::Disk*>(component->discretisation()), "Tried to access a hub with non-Disk discretisation!")

                        auto hub = dynamic_cast<discretisation::Disk*>(component->discretisation());
                        omega = hub->relativeRotationalVelocity()[2];

                    } // if blade

                    for( auto & child : children_ ) {
                        child.do_getHubOmega(omega);
                    }
                } // function do_getHubOmega

                HOST_PREFIX void do_getAABBInformation( real_t & towerHeight, Vector3<real_t> & hubPosition, real_t & hubRadius, real_t & bladeLength ) {

                    auto * discretisation = component->discretisation();
                    auto * points = discretisation->points();

                    if( !!(component->componentType() & ComponentType::TOWER) ) {
                        auto tower = dynamic_cast<discretisation::Line*>(discretisation);
                        WALBERLA_ASSERT_NOT_NULLPTR(tower)
                        towerHeight = (tower->endPoint().position - tower->startPoint().position).length();
                    } else if ( !!(component->componentType() & ComponentType::HUB) ) {
                        auto hub = dynamic_cast<discretisation::Disk*>(discretisation);
                        WALBERLA_ASSERT_NOT_NULLPTR(hub)
                        hubPosition = hub->points()[0].position;
                        hubRadius = hub->radius();
                    } else if ( !!(component->componentType() & ComponentType::BLADE) ) {
                        //TODO blade length is overwritten for every blade
                        bladeLength = (points[discretisation->nPoints()-1].position - points[0].position).length();
                    } else if ( !!(component->componentType() & ComponentType::ROTOR_DISK) ) {
                        auto rotorDisk = dynamic_cast<discretisation::Disk*>(discretisation);
                        WALBERLA_ASSERT_NOT_NULLPTR(rotorDisk)
                        bladeLength = rotorDisk->radius();
                    }

                    for( auto & child : children_ ) {
                        child.do_getAABBInformation( towerHeight, hubPosition, hubRadius, bladeLength );
                    }
                }

                /// MEMBER VARIABLES
                std::shared_ptr<Component> component{nullptr};

                Tree * parent_{nullptr};
                std::vector<Tree> children_{};

            };

        } // namespace tree

    } // namespace topology


} // namespace turbine_core

#endif //TURBINECORE_CPU_TREE_H
