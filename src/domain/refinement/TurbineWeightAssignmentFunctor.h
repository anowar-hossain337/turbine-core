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
//! \file TurbineWeightAssignmentFunctor.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef TURBINECORE_TURBINEWEIGHTASSIGNMENTFUNCTOR_H
#define TURBINECORE_TURBINEWEIGHTASSIGNMENTFUNCTOR_H

#include <blockforest/Block.h>
#include <blockforest/BlockForest.h>
#include <blockforest/loadbalancing/PODPhantomData.h>

#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace refinement {

        template< typename WindFarm_T >
        class TurbineWeightAssignmentFunctor {

        public:

            typedef walberla::blockforest::PODPhantomWeight<real_t> PhantomBlockWeight;

            explicit TurbineWeightAssignmentFunctor( WindFarm_T * const farm )
                    : farm_(farm)
            {}

            void operator()( std::vector< std::pair< const walberla::PhantomBlock *, walberla::any > > & blockData, const walberla::PhantomBlockForest & ){

                std::vector<walberla::AABB> turbineAABBs;
                farm_->calculateTurbineAABBs(turbineAABBs);

                for( auto &it : blockData ) {

                    const walberla::PhantomBlock * block = it.first;
                    WALBERLA_ASSERT_LESS( std::abs(int(block->getLevel()) - int(block->getSourceLevel())), 2 );

                    // basic LBM weight -> dependent on refinement level
                    real_t weight = real_t(1 << block->getLevel());

                    const walberla::AABB blockAABB = block->getAABB();

                    for(auto & aabb : turbineAABBs) {
                        if (blockAABB.intersects(aabb)) {
                            weight *= real_t(1.5);
                        }
                    }

                    it.second = PhantomBlockWeight( weight );

                }
            }

        private:

            WindFarm_T * const farm_;

        };

    } // namespace refinement

} // namespace turbine_core

#endif // TURBINECORE_TURBINEWEIGHTASSIGNMENTFUNCTOR_H