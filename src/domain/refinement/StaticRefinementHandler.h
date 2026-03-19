//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file TurbineAABBRefinementSelector.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_STATICREFINEMENTHANDLER_H
#define TURBINECORE_STATICREFINEMENTHANDLER_H

#pragma once

#include <utility>
#include <vector>

#include <core/config/Config.h>
#include <wind_turbine_core/math/AABB.h>

#include "wind_turbine_core/WalberlaDataTypes.h"

#include "domain/EnvironmentSetup.h"

namespace turbine_core {

    namespace refinement {

        enum BoundaryLayerPosition {
            NONE   = 1u << 0u,    // no boundary
            FRONT  = 1u << 1u,    // yMin
            BACK   = 1u << 2u,    // yMax
            EAST   = 1u << 3u,    // xMin
            WEST   = 1u << 4u,    // xMax
            TOP    = 1u << 5u,    // zMax
            BOTTOM = 1u << 6u     // zMin
        };

        template< typename WindFarm_T >
        class StaticRefinementHandler {

            using BLP_T = BoundaryLayerPosition;

        public:

            using AABBPair = std::pair< AABB, uint_t >;

            StaticRefinementHandler( const walberla::Config::BlockHandle & domainConfig, WindFarm_T * const farm,
                                     const EnvironmentSetup::Type & environmentSetup )
            : farm_(farm)
            {

                walberla::Config::BlockHandle refinementConfig = domainConfig;

                if( domainConfig.getNumBlocks("Refinement") ) {
                    refinementConfig = domainConfig.getOneBlock("Refinement");
                }

                finestLevel_ = refinementConfig.getParameter<uint_t>("numLevels", uint_t(0));

                if( refinementConfig.getNumBlocks("AABB") != 0 ) {

                    const auto & aabbConfig = refinementConfig.getOneBlock("AABB");

                    if( aabbConfig.getNumBlocks("turbine") != 0 ) {

                        std::vector<walberla::Config::BlockHandle> blocks{};
                        aabbConfig.getOneBlock("turbine").getBlocks(blocks);

                        for (const auto & block : blocks) {

                            const walberla::Vector3<real_t> min = block.getParameter<walberla::Vector3<real_t>>("min");
                            const walberla::Vector3<real_t> max = block.getParameter<walberla::Vector3<real_t>>("max");

                            const uint_t level = block.getParameter<uint_t>("level");
                            if( level > finestLevel_ ) {
                                WALBERLA_ABORT("Level of AABB must not be greater than finest level.")
                            }

                            turbineRelativeAABBs_.emplace_back(walberla::AABB(min,max), level);

                        }

                    }

                    if( aabbConfig.getNumBlocks("global") != 0 ) {

                        std::vector<walberla::Config::BlockHandle> blocks{};
                        aabbConfig.getOneBlock("global").getBlocks(blocks);

                        for (const auto & block : blocks) {

                            const walberla::Vector3<real_t> min = block.getParameter<walberla::Vector3<real_t>>("min");
                            const walberla::Vector3<real_t> max = block.getParameter<walberla::Vector3<real_t>>("max");

                            const uint_t level = block.getParameter<uint_t>("level");
                            if( level > finestLevel_ ) {
                                WALBERLA_ABORT("Level of AABB must not be greater than finest level.")
                            }

                            globalAABBs_.emplace_back(walberla::AABB(min,max), level);

                        }

                    }

                }

                // get initial turbine aabbs
                turbineAABBs_.clear();
                turbineAABBs_ = farm_->getTurbineAABBs();

                boundaryLayerLevel_ = refinementConfig.getParameter<uint_t>("BoundaryLayer", uint_t(0));

                // boundary layer
                if( environmentSetup == EnvironmentSetup::Open ) {
                    boundaryLayerPosition_ = BLP_T::BOTTOM;
                } else if ( environmentSetup == EnvironmentSetup::Tunnel ) {
                    boundaryLayerPosition_ = static_cast<BoundaryLayerPosition>(BLP_T::BOTTOM | BLP_T::TOP |
                                                                                BLP_T::BACK | BLP_T::FRONT);
                }

            } // constructor

            // for static mesh refinement
            void operator()( walberla::SetupBlockForest & forest ) {

                // calculate refinementAABBs
                std::vector< AABBPair > aabbs{};
                calculateRefinementAABBs(aabbs);

                if( aabbs.empty() )
                    return;

                for( auto block = forest.begin(); block != forest.end(); ++block ) {

                    if( boundaryLayerPosition_ != BLP_T::NONE && boundaryLayerLevel_ != 0 ) {
                        // add boundary layer
                        if (checkIfBoundary(block.get(),forest) && block->getLevel() < boundaryLayerLevel_) {
                            block->setMarker(true);
                        }
                    }

                    // add aabbs
                    for(auto & aabb : aabbs) {
                        if( aabb.first.intersects(block->getAABB()) && block->getLevel() < aabb.second )
                            block->setMarker( true );
                    }
                }
            }

            void calculateRefinementAABBs( std::vector<AABBPair> & refinementAABBs ) {

                refinementAABBs.clear();

                // add global AABBs
                for( const auto & aabb : globalAABBs_ ) {
                    refinementAABBs.emplace_back(aabb);
                }

                turbineAABBs_.clear();
                turbineAABBs_ = farm_->getTurbineAABBs();

                // calculate and add turbine relative AABBs
                for( const auto & turbineAABB : turbineAABBs_ ) {

                    refinementAABBs.emplace_back(turbineAABB, finestLevel_);

                    const auto hubPosition = turbineAABB.center();

                    for( const auto & relAABB : turbineRelativeAABBs_ ) {

                        auto centeredAABB = AABB( hubPosition + relAABB.first.min(),
                                                  hubPosition + relAABB.first.max());

                        refinementAABBs.emplace_back( centeredAABB, relAABB.second );
                    }

                }

            }

            // for dynamic refinement
            void operator()(std::vector< std::pair< const walberla::Block *, uint_t > > & minTargetLevels,
                            std::vector< const walberla::Block * > &, const walberla::BlockForest & forest ) {

                // calculate refinementAABBs
                std::vector< AABBPair > aabbs{};
                calculateRefinementAABBs(aabbs);

                for( auto it = minTargetLevels.begin(); it != minTargetLevels.end(); ++it ) {
                    uint_t currentLevelOfBlock = it->first->getLevel();
                    uint_t targetLevelOfBlock = currentLevelOfBlock;

                    if( boundaryLayerPosition_ != BLP_T::NONE && boundaryLayerLevel_ != 0 ) {
                        // add boundary layer
                        if (checkIfBoundary(it->first,forest)) {
                            // set block to finest level since overlap with at least one global body is present
                            uint_t targetLevelOfBoundary = boundaryLayerLevel_;

                            if (currentLevelOfBlock > targetLevelOfBoundary) {
                                targetLevelOfBlock = currentLevelOfBlock - uint_t(1);
                            } else if (currentLevelOfBlock < targetLevelOfBlock) {
                                targetLevelOfBlock = currentLevelOfBlock + uint_t(1);
                            }

                            // only the first found intersecting AABB is taken into account
                            break;
                        }
                    }

                    // add aabbs
                    for(auto & aabb : aabbs) {

                        if( aabb.first.intersects(it->first->getAABB()) ) {

                            uint_t targetLevelOfAABB = aabb.second;

                            if( currentLevelOfBlock > targetLevelOfAABB ) {
                                targetLevelOfBlock = currentLevelOfBlock - uint_t(1);
                            } else if ( currentLevelOfBlock < targetLevelOfBlock ) {
                                targetLevelOfBlock = currentLevelOfBlock + uint_t(1);
                            }

                            // only the first found intersecting AABB is taken into account
                            break;

                        }
                    }

                    WALBERLA_CHECK_LESS_EQUAL(std::abs(int(targetLevelOfBlock) - int(currentLevelOfBlock)), uint_t(1), "Only level difference of maximum 1 allowed!");
                    it->second = targetLevelOfBlock;
                }

            }

            WindFarm_T * getFarm() const {
                return farm_;
            }

            auto getTurbineAABBs() const {
                return turbineAABBs_;
            }

            template<typename Block_T, typename Forest_T>
            bool checkIfBoundary(const Block_T * block, Forest_T & bf) {

                return (boundaryLayerPosition_ & BLP_T::EAST   && bf.atDomainXMinBorder(*block)) ||
                       (boundaryLayerPosition_ & BLP_T::WEST   && bf.atDomainXMaxBorder(*block)) ||
                       (boundaryLayerPosition_ & BLP_T::FRONT  && bf.atDomainYMinBorder(*block)) ||
                       (boundaryLayerPosition_ & BLP_T::BACK   && bf.atDomainYMaxBorder(*block)) ||
                       (boundaryLayerPosition_ & BLP_T::BOTTOM && bf.atDomainZMinBorder(*block)) ||
                       (boundaryLayerPosition_ & BLP_T::TOP    && bf.atDomainZMaxBorder(*block));

            }

        private:

            WindFarm_T * const farm_{};

            std::vector<AABB> turbineAABBs_{};
            std::vector<AABBPair> turbineRelativeAABBs_{};
            std::vector<AABBPair> globalAABBs_{};
            uint_t finestLevel_{};

            BoundaryLayerPosition boundaryLayerPosition_{BoundaryLayerPosition::NONE};
            uint_t boundaryLayerLevel_{};

        };

    } // namespace refinement

} // namespace turbine_core

#endif // TURBINECORE_STATICREFINEMENTHANDLER_H
