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
//! \file TurbineSelectorResetter.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef TURBINECORE_TURBINESELECTORRESETTER_H
#define TURBINECORE_TURBINESELECTORRESETTER_H

#include <blockforest/BlockForest.h>

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/Selectors.h"

namespace turbine_core {

    namespace refinement {

        template< typename WindTurbine_T >
        class TurbineSelectorResetter {

        public:

            TurbineSelectorResetter( WindTurbine_T * const farm, const uint_t & ghostLayers )
                    : farm_(farm), ghostLayers_(ghostLayers)
            {}

            void operator()( walberla::BlockForest & forest, const walberla::PhantomBlockForest & phantomForest );

        private:

            WindTurbine_T * const farm_;
            const uint_t ghostLayers_;

        };

        template< typename WindTurbine_T >
        void TurbineSelectorResetter<WindTurbine_T>::operator()(walberla::BlockForest & forest, const walberla::PhantomBlockForest & ) {

            std::vector<walberla::AABB> turbineAABBs;
            farm_->calculateTurbineAABBs(turbineAABBs);

            // clear boundary handling & reset turbine flags
            for (auto block = forest.begin(); block != forest.end(); ++block) {

                // clear block flags
                block->clearState();

                // set new state - turbines
                for (uint_t t = 0; t < turbineAABBs.size(); ++t) {

                    auto blockAABB = block->getAABB().getExtended(ghostLayers_);

                    if (blockAABB.intersects(turbineAABBs[t])) {
                        block->addState(sets::turbineSelector());
                        block->addState(sets::turbineSelector(t));
                    }
                }

                // set new state - refinement
                if (forest.getLevel(*block) > 0) {
                    block->addState(sets::refinementSelector());
                }

            }

        }

    } // namespace refinement

} // namespace turbine_core

#endif // TURBINECORE_TURBINESELECTORRESETTER_H