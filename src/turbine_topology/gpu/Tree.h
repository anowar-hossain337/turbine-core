
#pragma once

#ifndef TURBINECORE_GPU_TREE_H
#define TURBINECORE_GPU_TREE_H

#include "turbine_topology/gpu/GPUTurbineTopology.h"
#include "turbine_topology/Stack.h"
#include "component/Component.h"
#include "conversion/Conversion.h"
#include "wind_turbine_core/ProjectDefines.h"

#include <cmath>
#include <cassert>

namespace turbine_core {

    namespace topology {

        namespace gpu {

            class Tree : public GPUTurbineTopology<Tree> {

                friend GPUTurbineTopology<Tree>;

            public:

                HOST_DEVICE_PREFIX Tree() {}

                HOST_DEVICE_PREFIX Tree( const char * ID, const ComponentType & type, Tree * parent,
                                         discretisation::Discretisation * discretisation,
                                         force_model::ForceModel * forceModel = nullptr,
                                         control_model::ControlModel * controlModel = nullptr,
                                         Point3<real_t> * reference = nullptr )
                        : parent_(parent)
                {

                    createListFromTree();
                    assert(treeList_ != nullptr);

                    //TODO do I want this or an assert?
                    if(discretisation != nullptr) {
                        delete component;
                        component = new Component(ID, type, discretisation, forceModel, controlModel, reference);
                    }

                }

                HOST_DEVICE_PREFIX Tree( const Tree & tree )
                : parent_(tree.parent_), nChildren_(tree.nChildren_), capacity_(tree.capacity_) {

                    //TODO meehhh....
                    if(tree.component != nullptr) {
                        delete component;
                        component = new Component(*tree.component);
                    }

                    free(children_);
                    if(tree.capacity_ != 0) {
			            children_ = (Tree*) malloc( tree.capacity_ * sizeof(Tree) );
		            } else {
                        children_ = nullptr;
                    }

                    for(size_t i = 0; i < tree.nChildren_; ++i) {
                        new (static_cast<Tree*>(children_) + i ) Tree(tree.children_[i]);
                    }

                    createListFromTree();

                }

                HOST_DEVICE_PREFIX Tree& operator=( const Tree & tree ) {
                    parent_ = tree.parent_;
                    nChildren_ = tree.nChildren_;
                    capacity_ = tree.capacity_;

                    //TODO meehhh....
                    if(tree.component != nullptr) {
                        delete component;
                        component = new Component(*tree.component);
                    }

                    free(children_);
                    if(tree.capacity_ != 0) {
                        children_ = (Tree*) malloc( tree.capacity_ * sizeof(Tree) );
                    } else {
                        children_ = nullptr;
                    }

                    for(size_t i = 0; i < tree.nChildren_; ++i) {
                        new (static_cast<Tree*>(children_) + i ) Tree(tree.children_[i]);
                    }

                    createListFromTree();

                    return *this;
                }

                HOST_DEVICE_PREFIX ~Tree() {

                    free(treeList_);
                    delete component;

                    for(uint_t i = 0; i < nChildren_; ++i) {
                        (static_cast<Tree*>(children_) + i )->~Tree();
                    }
                    free(children_);
                }

                HOST_DEVICE_PREFIX static bool compareID( const char * ID1, const char * ID2 ) {

                    auto * p1 = (const unsigned char *) ID1;
                    auto * p2 = (const unsigned char *) ID2;

                    while (*p1 && (*p1 == *p2)) {
                        ++p1, ++p2;
                    }

                    return ((*p1 > *p2) - (*p2 > *p1)) == 0;
                }

                HOST_DEVICE_PREFIX void createListFromTree() {

                    free(treeList_);
                    treeList_ = (Tree**)malloc(15*sizeof(Tree*)); // NOLINT(bugprone-sizeof-expression)
                    auto * s = stack::createStack<Tree>(15);
                    push(s, this);
                    int nComponents = 0;
                    Tree * temp;
                    while(!isEmpty(s)){
                        temp = pop(s);
                        treeList_[nComponents++] = temp;

                        for( size_t c = temp->nChildren_; c > 0; --c ) {
                            push(s, &temp->children_[c-1]);
                        }
                    }

                    treeLength_ = nComponents;

                    free(s->node);
                    free(s);
                }

                HOST_DEVICE_PREFIX void do_callback( const Component::Function & function, const ComponentType & type, const uint_t timestep) const {
                    for(int c = 1; c < treeLength_; ++c) {
                        treeList_[c]->component->callFunction(function, type, timestep);
                    }
                }

                HOST_DEVICE_PREFIX void do_callback( const Component::Output & function, const ComponentType & type ) const {
                    for(int c = 1; c < treeLength_; ++c) {
                        treeList_[c]->component->callFunction(function, type);
                    }
                }

                HOST_DEVICE_PREFIX void do_addDiscretisation (const char * parentID, const char * ID,
                                                              const ComponentType & type,
                                                              discretisation::Discretisation * discretisation ) {

                    for(int c = 0; c < treeLength_; ++c) {

                        Tree * current = treeList_[c];

                        if (compareID(parentID, current->component->componentID())) {

                            if (current->nChildren_ == current->capacity_) {

                                current->capacity_ = max(current->capacity_ + 1, 2 * current->capacity_);

                                Tree * tmp = (Tree *) malloc(current->capacity_ * sizeof(Tree));

                                for (size_t i = 0; i < current->nChildren_; ++i) {
                                    new (static_cast<Tree*>(tmp) + i ) Tree(current->children_[i]);
                                    (current->children_[i]).~Tree();
                                }

                                free(current->children_);
                                current->children_ = tmp;

                            }

                            auto * reference = current->component->getReferencePointForChild();

                            new (static_cast<Tree*>(current->children_) + current->nChildren_ ) Tree(ID, type, current, discretisation, nullptr, nullptr, reference);

                            ++(current->nChildren_);

                            createListFromTree();

                            return;
                        }
                    }

                    printf("Please insert valid parent ID.\n");

                }

                HOST_DEVICE_PREFIX void do_initialiseForceAndControlModel (const char * ID, force_model::ForceModel * forceModel, control_model::ControlModel * controlModel ) {

                    for(int c = 0; c < treeLength_; ++c) {

                        Tree * current = treeList_[c];

                        if (compareID(ID, current->component->componentID())) {
                            current->component->initialiseForceAndControlModel(forceModel, controlModel);
                            return;
                        }
                    }

                    printf("No component with this ID found when trying to add force model.\n");

                }

                HOST_DEVICE_PREFIX void do_print() const {
                    for(int c = 1; c < treeLength_; ++c) {
                        printf("\nComponent: %s (type = %iu)\n", treeList_[c]->component->componentID(), static_cast<unsigned int>(treeList_[c]->component->componentType()));
                        treeList_[c]->component->print();
                    }
                }

                HOST_DEVICE_PREFIX void do_applyControl(const uint_t & timeStep, const real_t & rotorMeanWind, const real_t & rotorTorque,
                    const real_t rotorOmega, const Vector3<real_t> rotorWindVaneWind) {
                    for(int c = 1; c < treeLength_; ++c) {
                        treeList_[c]->component->applyControl(timeStep, rotorMeanWind, rotorTorque, rotorOmega, rotorWindVaneWind);
                    }
                }

                HOST_DEVICE_PREFIX void do_getControlNeeds(bool & needsMeanVelocity, bool & needsTorque, bool & needsWindVane)
                {
                    for(int c = 1; c < treeLength_; ++c) {
						bool meanVel{false}, torque{false}, windVane{false};
                        treeList_[c]->component->getControlNeeds(meanVel, torque, windVane);
						needsMeanVelocity += meanVel;
						needsTorque += torque;
						needsWindVane += windVane;
                    }
                }

                HOST_DEVICE_PREFIX void do_evaluateDensityAndVelocity( const blockforest::BlockInfo & blockInfo,
                                                                       const projectors::Projectors &projectors,
                                                                       uint_t pointIdx, uint_t pointIdy)
                {
                    if (pointIdx < treeLength_) {
                        treeList_[pointIdx]->component->evaluateDensityAndVelocity(blockInfo, projectors, pointIdy);
                    }
                }

                HOST_DEVICE_PREFIX void do_spreadForces( const blockforest::BlockInfo & blockInfo,
                                                         const projectors::Projectors &projectors,
                                                         uint_t pointIdy)
                {
                    for(int c = 1; c < treeLength_; ++c) {
                        treeList_[c]->component->spreadForces(blockInfo, projectors, pointIdy);
                    }
                }

                HOST_DEVICE_PREFIX void do_setupCommunication(const uint_t ownRank, const uint_t neighborRank,
                                                              int * sendPoints, int * recvPoints,
                                                              const uint_t nLocalAABBs, const uint_t nNeighborAABBs,
                                                              math::AABB * const localAABBs, math::AABB * const neighborAABBs,
                                                              uint_t & nSendPoints, walberla::mpi::MPISize & receiverSize,
                                                              const mpi::SynchronisationMode mode, const real_t impactWidth) {

                    // reset variables in case they contain stupid stuff
                    nSendPoints = uint_t(0);
                    receiverSize = walberla::mpi::MPISize(0);

                    uint_t nBlade{0};
                    const real_t impactWidthSq = real_t(4.0) * impactWidth * impactWidth;

                    for(int c = 1; c < treeLength_; ++c) {

                        if (!!(treeList_[c]->component->componentType() & ComponentType::BLADE)) {

                            auto * discretisation = treeList_[c]->component->discretisation();
                            assert(discretisation && "Tried to access a blade with non-Line discretisation!");

                            const uint_t nPoints = discretisation->nPoints();
                            const auto *points = discretisation->points();

                            int * send = &sendPoints[nBlade * (nPoints+1)];
                            int * recv = &recvPoints[nBlade * (nPoints+1)];

                            uint_t nSend{0};
                            uint_t nRecv{0};

                            // loop over all points to get data size to communicate
                            for (uint_t p = 0; p < nPoints; ++p) {
                                const auto pos = points[p].position;

                                bool containedLocal{false};
                                bool intersectsLocal{false};
                                for(uint_t n = 0 ; n < nLocalAABBs; ++n) {
                                    if(domain::sqDistancePointToAABB(pos, localAABBs[n]) <= real_t(0)) {
                                        containedLocal = true;
                                    } else if(domain::sqDistancePointToAABB(pos, localAABBs[n]) <= impactWidthSq) {
                                        intersectsLocal = true;
                                    }
                                    if(containedLocal && intersectsLocal)
                                        break;
                                }

                                bool containedRemote{false};
                                bool intersectsRemote{false};
                                for(uint_t n = 0 ; n < nNeighborAABBs; ++n) {
                                    if(domain::sqDistancePointToAABB(pos, neighborAABBs[n]) <= real_t(0)) {
                                        containedRemote = true;
                                    } else if(domain::sqDistancePointToAABB(pos, neighborAABBs[n]) <= impactWidthSq) {
                                        intersectsRemote = true;
                                    }
                                    if(containedRemote && intersectsRemote)
                                        break;
                                }

                                // point lies in local block but impact width reaches into remote block
                                //   -> send to remote neighbour
                                if (containedLocal && intersectsRemote) {
                                    assert(p < internal::nBladePoints);
                                    send[nSend] = int(p);
                                    ++nSend;
                                }

                                // point lies in remote block but impact width reaches into local block
                                //   -> receive from neighbour
                                if (intersectsLocal && containedRemote) {
                                    assert(p < internal::nBladePoints);
                                    recv[nRecv] = int(p);
                                    ++nRecv;
                                    receiverSize += walberla::mpi::MPISize(mpi::itemSize(mode));
                                }
                            }

                            // assure that other points are -1
                            for( uint_t p = nSend; p < nPoints+1; ++p ) {
                                send[p] = -1;
                            }
                            for( uint_t p = nRecv; p < nPoints+1; ++p ) {
                                recv[p] = -1;
                            }

                            nSendPoints += nSend;
                            ++nBlade;
                        }

                        if (!!(treeList_[c]->component->componentType() & ComponentType::ROTOR_DISK)) {

                            auto * discretisation = treeList_[c]->component->discretisation();
                            assert(discretisation && "Tried to access a blade with non-Line discretisation!");

                            const uint_t nPoints = discretisation->nPoints();
                            const auto *points = discretisation->points();

                            int * send = &sendPoints[0];
                            int * recv = &recvPoints[0];

                            uint_t nSend{0};
                            uint_t nRecv{0};

                            // loop over all points to get data size to communicate
                            for (uint_t p = 0; p < nPoints; ++p) {
                                const auto pos = points[p].position;

                                bool containedLocal{false};
                                bool intersectsLocal{false};
                                for(uint_t n = 0 ; n < nLocalAABBs; ++n) {

                                    auto sqDistance = domain::sqDistancePointToAABB(pos, localAABBs[n]);
                                    if(sqDistance <= real_t(0)) {
                                        containedLocal = true;
                                    } else if(sqDistance <= impactWidthSq) {
                                        intersectsLocal = true;
                                    }
                                }

                                bool containedRemote{false};
                                bool intersectsRemote{false};
                                for(uint_t n = 0 ; n < nNeighborAABBs; ++n) {

                                    auto sqDistance = domain::sqDistancePointToAABB(pos, neighborAABBs[n]);
                                    if(sqDistance <= real_t(0)) {
                                        containedRemote = true;
                                    } else if(sqDistance  <= impactWidthSq) {
                                        intersectsRemote = true;
                                    }
                                }
                                // point lies in local block but impact width reaches into remote block
                                //   -> send to remote neighbour
                                if (containedLocal && intersectsRemote) {
                                    send[nSend] = int(p);
                                    ++nSend;
                                }

                                // point lies in remote block but impact width reaches into local block
                                //   -> receive from neighbour
                                if (intersectsLocal && containedRemote) {
                                    recv[nRecv] = int(p);
                                    ++nRecv;
                                    receiverSize += walberla::mpi::MPISize(mpi::itemSize(mode));
                                }
                            }

                            // assure that other points are -1
                            for( uint_t p = nSend; p < nPoints+1; ++p ) {
                                send[p] = -1;
                            }
                            for( uint_t p = nRecv; p < nPoints+1; ++p ) {
                                recv[p] = -1;
                            }

                            nSendPoints += nSend;
                        }

                    }

                }

                template<typename Type_T>
                HOST_DEVICE_PREFIX void do_packTurbineData(Type_T * buffer,
                                                           int const * const sendPoints,
                                                           const mpi::SynchronisationMode mode) {
                    uint_t nBlade{0};
                    for (int c = 1; c < treeLength_; ++c) {

                        if (!!(treeList_[c]->component->componentType() & ComponentType::BLADE)) {

                            auto * idx = &sendPoints[nBlade * (internal::nBladePoints+1)];

                            while( *idx != -1 ) {

                                treeList_[c]->component->packData(buffer, mode, *idx);

                                buffer += static_cast<size_t>(mpi::itemSize(mode));
                                ++idx;
                            }

                            ++nBlade;
                        }

                        if (!!(treeList_[c]->component->componentType() & ComponentType::ROTOR_DISK)) {
                            auto * idx = &sendPoints[0];

                            while( *idx != -1 ) {

                                treeList_[c]->component->packData(buffer, mode, *idx);

                                buffer += static_cast<size_t>(mpi::itemSize(mode));
                                ++idx;
                            }
                        }
                    }
                }

                template<typename Type_T>
                HOST_DEVICE_PREFIX void do_unpackTurbineData(Type_T * buffer,
                                                             int * recvPoints,
                                                             const mpi::SynchronisationMode mode) {

                    uint_t nBlade{0};
                    for (int c = 1; c < treeLength_; ++c) {

                        if (!!(treeList_[c]->component->componentType() & ComponentType::BLADE)) {

                            auto * idx = &recvPoints[nBlade * (internal::nBladePoints+1)];

                            // storing only uint_t's for indices, rest is initialised with -1
                            while( *idx != -1 ) {

                                treeList_[c]->component->unpackData(buffer, mode, *idx);

                                buffer += static_cast<size_t>(mpi::itemSize(mode));
                                ++idx;
                            }

                            ++nBlade;
                        }

                        if (!!(treeList_[c]->component->componentType() & ComponentType::ROTOR_DISK)) {

                            auto * idx = &recvPoints[0];

                            // storing only uint_t's for indices, rest is initialised with -1
                            while( *idx != -1 ) {

                                treeList_[c]->component->unpackData(buffer, mode, *idx);

                                buffer += static_cast<size_t>(mpi::itemSize(mode));
                                ++idx;
                            }
                        }
                    }
                }

                HOST_DEVICE_PREFIX void do_getBladeData() {
                    internal::nBlades = 0;
                    for(int c = 1; c < treeLength_; ++c) {
                        if(!!(treeList_[c]->component->componentType() & ComponentType::BLADE)) {
                            ++internal::nBlades;
                            internal::nBladePoints = treeList_[c]->component->discretisation()->nPoints();
                        }
                    }
                }

                HOST_DEVICE_PREFIX uint_t do_getNComponents() const {
                    return (uint_t)treeLength_;
                }

                HOST_DEVICE_PREFIX uint_t do_getNControlledComponents() const {
                    uint_t nControlledComponents;
                    for( uint_t c = 0; c < treeLength_; ++c ) {
                        if(treeList_[c]->component->isControlled()){
                            nControlledComponents++;
                        }
                    }
                    return nControlledComponents;
                }

                HOST_DEVICE_PREFIX void do_getNPoints( uint_t * nPoints ) {
                    for( uint_t c = 0; c < treeLength_; ++c ) {
                        nPoints[c] = treeList_[c]->component->discretisation()->nPoints();
                    }
                }

                template< typename DataType_T >
                HOST_DEVICE_PREFIX void do_getOutputData( const Component::Output & function, DataType_T * data ) const {

                    uint_t nDataItems{1};
                    if(function == Component::Output::ORIENTATIONS)
                        nDataItems *= 4;

                    uint_t currentIdx = 0;
                    for(int c = 0; c < treeLength_; ++c) {

                        auto discretisation = treeList_[c]->component->discretisation();
                        auto points = discretisation->points();

                        for(uint_t p = 0; p < discretisation->nPoints(); ++p) {

                            const uint_t idx = (currentIdx + p) * nDataItems;

                            data[idx] = (points+p)->position;

                            if(function == Component::Output::ORIENTATIONS) {
                                auto matrix = (points+p)->orientation.toRotationMatrix();
                                data[idx+1] = Vector3<real_t>(matrix[0], matrix[3], matrix[6]);
                                data[idx+2] = Vector3<real_t>(matrix[1], matrix[4], matrix[7]);
                                data[idx+3] = Vector3<real_t>(matrix[2], matrix[5], matrix[8]);
                            }

                        }

                        currentIdx += discretisation->nPoints();
                    }
                }

                HOST_DEVICE_PREFIX void do_getForceData( const blockforest::BlockInfo & blockInfo,
                                                         uint_t * const nLocalPoints, uint_t * const localPoints,
                                                         real_t * const azimuth, Vector2<real_t> * const localForces,
                                                         Vector3<real_t> * const localAero, real_t * const relVelocity) {

                    uint_t currentBlade{0};
                    *nLocalPoints = 0;

                    for(int c = 1; c < treeLength_; ++c) {
                        if(!!(treeList_[c]->component->componentType() & ComponentType::BLADE)) {

                            //TODO might be error prone... what if discretisation is sth else?
                            auto * line = static_cast<discretisation::Line*>(treeList_[c]->component->discretisation());
                            assert(line && "Tried to access a blade with non-Line discretisation!");

                            const uint_t nPoints = line->nPoints();

                            auto reference = treeList_[c]->parent_->component->getReferencePointFromParent();

                            auto * discretisationPoints = treeList_[c]->component->discretisation()->points();
                            auto * forcePoints = treeList_[c]->component->forceModel()->points();

                            assert(discretisationPoints && "Component returned null pointer as discretisation points.");
                            assert(forcePoints && "Component returned null pointer as actuator data points.");

                            const real_t factor = real_t(1.0) / line->elementLength_;

                            const auto & aabb = blockInfo.getAABB();

                            for(uint_t p = 0; p < line->nPoints(); ++p) {
                                if(!aabb.contains(discretisationPoints[p].position))
                                    continue;

                                line->calculateAzimuthAngle(reference, p, azimuth[*nLocalPoints]);

                                Vector3<real_t> localForces3 = discretisationPoints[p].orientation.getInverse().rotate(forcePoints[p].forcesGlobal);
                                localForces[*nLocalPoints] = Vector2<real_t>(localForces3[0] * factor, localForces3[1] * factor);

                                // Additional variables
                                auto aoa = forcePoints[p].angleOfAttack;
                                auto polars = forcePoints[p].polar(aoa);

                                localAero[*nLocalPoints] = Vector3<real_t>(aoa, polars[0], polars[1]);

                                auto relativeVelocity = forcePoints[p].velocityLocal;

                                relativeVelocity[2] = 0.;
                                relVelocity[*nLocalPoints] = relativeVelocity.length();

                                localPoints[*nLocalPoints] = currentBlade * nPoints + p;
                                (*nLocalPoints) += 1;
                            }

                            ++currentBlade;
                        } // if blade
                    } // loop components
                }

                HOST_DEVICE_PREFIX void do_calculatePowerAndThrust( const blockforest::BlockInfo & blockInfo, real_t & power, real_t & thrust ) {

                    power = thrust = real_t(0);

                    for(int c = 1; c < treeLength_; ++c) {
                        if(!!(treeList_[c]->component->componentType() & ComponentType::BLADE)) {

                            auto * line = static_cast<discretisation::Line*>(treeList_[c]->component->discretisation());
                            assert(line && "Tried to access a blade with non-Line discretisation!");

                            auto * discretisationPoints = treeList_[c]->component->discretisation()->points();
                            assert(discretisationPoints && "Component returned null pointer as discretisation points.");

                            auto * forcePoints = treeList_[c]->component->forceModel()->points();
                            assert(forcePoints && "Component returned null pointer as actuator data points.");

                            const uint_t nPoints = treeList_[c]->component->discretisation()->nPoints();
                            auto reference = treeList_[c]->component->getReferencePointFromParent();

                            auto orientation = reference->orientation;

                            auto hub = static_cast<discretisation::Disk*>(treeList_[c]->parent_->component->discretisation());
                            auto rotVelocity = hub->relativeRotationalVelocity();

                            const auto& blockBB = blockInfo.getAABB();

                            for(uint_t p = 0; p < nPoints; ++p) {
                                const auto position = discretisationPoints[p].position;
                                if( blockBB.contains(position) ) {
                                    // thrust and power are local
                                    thrust += orientation.getInverse().rotate(forcePoints[p].forcesGlobal)[2];
                                    power += orientation.getInverse().rotate((position - reference->position) % forcePoints[p].forcesGlobal)[2] * rotVelocity[2];
                                }
                            }
                        } // if blade
                    } // loop components
                }

                HOST_DEVICE_PREFIX void do_calculateTorque( const blockforest::BlockInfo & blockInfo, real_t & torque ) {

                    torque = real_t(0);

                    for(int c = 1; c < treeLength_; ++c) {
                        if(!!(treeList_[c]->component->componentType() & ComponentType::BLADE)) {

                            auto * line = static_cast<discretisation::Line*>(treeList_[c]->component->discretisation());
                            assert(line && "Tried to access a blade with non-Line discretisation!");

                            auto * discretisationPoints = treeList_[c]->component->discretisation()->points();
                            assert(discretisationPoints && "Component returned null pointer as discretisation points.");

                            auto * forcePoints = treeList_[c]->component->forceModel()->points();
                            assert(forcePoints && "Component returned null pointer as actuator data points.");

                            const uint_t nPoints = treeList_[c]->component->discretisation()->nPoints();
                            auto reference = treeList_[c]->component->getReferencePointFromParent();

                            auto orientation = reference->orientation;

                            auto hub = static_cast<discretisation::Disk*>(treeList_[c]->parent_->component->discretisation());
                            auto rotVelocity = hub->relativeRotationalVelocity();

                            const auto& blockBB = blockInfo.getAABB();

                            for(uint_t p = 0; p < nPoints; ++p) {
                                const auto position = discretisationPoints[p].position;
                                if( blockBB.contains(position) ) {
                                    // torque is local
                                    torque += orientation.getInverse().rotate((position - reference->position) % forcePoints[p].forcesGlobal)[2];
                                }
                            }
                        } // if blade
                    } // loop components
                }

                HOST_DEVICE_PREFIX void do_getOmega( real_t & omega ) {

                    omega = real_t(0);

                    for(int c = 1; c < treeLength_; ++c) {
                        if(!!(treeList_[c]->component->componentType() & ComponentType::HUB)) {

                            auto * hub = static_cast<discretisation::Disk*>(treeList_[c]->component->discretisation());
                            assert(hub && "Tried to access a hub with non-Disk discretisation!");
                            omega = hub->relativeRotationalVelocity()[2];
                        } // if blade
                    } // loop components
                }

                HOST_DEVICE_PREFIX void do_calculateMeanVelocity( const blockforest::BlockInfo & blockInfo, real_t & totalVelocity, uint_t & nbPoints ) {

                    totalVelocity = real_t(0);
                    nbPoints = uint_t(0);

                    for(int c = 1; c < treeLength_; ++c) {
                        if(!!(treeList_[c]->component->componentType() & ComponentType::BLADE)) {

                            auto * discretisationPoints = treeList_[c]->component->discretisation()->points();
                            assert(discretisationPoints && "Component returned null pointer as discretisation points.");

                            auto * forceModel = treeList_[c]->component->forceModel();
                            assert(forceModel && "No force model assigned to blades.");

                            auto * forcePoints = forceModel->points();
                            assert(forcePoints && "Component returned null pointer as actuator data points.");

                            const uint_t nPoints = treeList_[c]->component->discretisation()->nPoints();
                            auto reference = treeList_[c]->component->getReferencePointFromParent();

                            auto orientation = reference->orientation;

                            const auto& blockBB = blockInfo.getAABB();

                            for(uint_t p = 0; p < nPoints; ++p) {
                                const auto position = discretisationPoints[p].position;
                                if( blockBB.contains(position) ) {
                                    ++nbPoints;
                                    // totalVelocity is local
                                    totalVelocity -= orientation.getInverse().rotate(forcePoints[p].velocityLocal)[2] ; // Velocity along the hub axis
                                }
                            }
                        } // if blade
                    } // loop components
                }

                HOST_DEVICE_PREFIX void do_calculateWindVaneVelocity( const blockforest::BlockInfo & blockInfo, Vector3<real_t> & totalVelocity, uint_t & nbPoints ) {

                    totalVelocity = Vector3<real_t> (0,0,0);

                    for(int c = 1; c < treeLength_; ++c) {
                        if(!!(treeList_[c]->component->componentType() & ComponentType::BLADE)) {

                            auto * discretisationPoints = treeList_[c]->component->discretisation()->points();
                            assert(discretisationPoints && "Component returned null pointer as discretisation points.");

                            auto * forceModel = treeList_[c]->component->forceModel();
                            assert(forceModel && "No force model assigned to blades.");

                            auto * forcePoints = forceModel->points();
                            assert(forcePoints && "Component returned null pointer as actuator data points.");

                            const uint_t nPoints = treeList_[c]->component->discretisation()->nPoints();
                            auto reference = treeList_[c]->component->getReferencePointFromParent();

                            auto orientation = reference->orientation;

                            const auto& blockBB = blockInfo.getAABB();

                            const auto position = discretisationPoints[0].position;
                            if( blockBB.contains(position)) {
                                ++nbPoints;
                                // totalVelocity is global
                                totalVelocity += forcePoints[0].velocityLocal; // Velocity along the hub axis
                            }
                        } // if blade
                    } // loop components
                }

                HOST_DEVICE_PREFIX void do_getAABBInformation( real_t * info ) {

                    real_t bladeLength, hubRadius, towerHeight;
                    Vector3<real_t> hubPosition{};

                    for(int c = 1; c < treeLength_; ++c) {

                        if(!!(treeList_[c]->component->componentType() & ComponentType::TOWER)) {

                            auto * tower = static_cast<discretisation::Line*>(treeList_[c]->component->discretisation());
                            assert(tower && "Tried to access a tower with non-Line discretisation!");

                            towerHeight = (tower->endPoint().position - tower->startPoint().position).length();

                        } else if(!!(treeList_[c]->component->componentType() & ComponentType::HUB)) {

                            auto hub = static_cast<discretisation::Disk*>(treeList_[c]->component->discretisation());
                            assert(hub && "Tried to access a hub with non-disk discretisation!");
                            hubPosition = hub->points()[0].position;
                            hubRadius = hub->radius();

                        } else if(!!(treeList_[c]->component->componentType() & ComponentType::BLADE)) {
                            auto * blade = static_cast<discretisation::Line*>(treeList_[c]->component->discretisation());
                            bladeLength = (blade->endPoint().position - blade->startPoint().position).length();
                        } else if(!!(treeList_[c]->component->componentType() & ComponentType::ROTOR_DISK)) {
                            auto rotorDisk = static_cast<discretisation::Disk*>(treeList_[c]->component->discretisation());				
                            bladeLength = rotorDisk->radius();
			            }
                    } // loop components

                    info[0] = bladeLength;
                    info[1] = hubRadius;
                    info[2] = towerHeight;
                    info[3] = hubPosition[0];
                    info[4] = hubPosition[1];
                    info[5] = hubPosition[2];

                }

            private:

                /// MEMBER VARIABLES
                component::Component * component{nullptr};

                Tree * parent_{nullptr};

                size_t nChildren_{0};
                size_t capacity_{0};
                Tree * children_{nullptr};

                int treeLength_{0};
                Tree ** treeList_{nullptr};

            };

        } // namespace gpu

    } // namespace data_structure

} // namespace turbine_core

#endif //TURBINECORE_GPU_TREE_H
