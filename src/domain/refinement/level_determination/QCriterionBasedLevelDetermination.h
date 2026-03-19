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
//! \file QCriterionBasedLevelDetermination.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_QCRITERIONBASEDLEVELDETERMINATION_H
#define TURBINECORE_QCRITERIONBASEDLEVELDETERMINATION_H

#pragma once

#include <blockforest/Block.h>
#include <core/math/Vector3.h>
#include <lbm/field/QCriterion.h>

#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace refinement {

        /*
         * Refinement check based on Q criterion
         * If Q criterion is below lowerLimit in all cells of a block, that block could be coarsened.
         * If the Q criterion is above the upperLimit for at least one cell, that block gets marked for refinement.
         * Else, the block remains on the current level.
         */
        template< typename Stencil_T, typename VectorField_T, typename Filter_T >
        class QCriterionBasedLevelDetermination
        {
        public:

            QCriterionBasedLevelDetermination( const BlockDataID & fieldID, const Filter_T & filter, walberla::StructuredBlockStorage * const storage,
                                               const real_t upperLimit, const real_t lowerLimit, const uint_t maxLevel ) :
                    fieldID_( fieldID ), filter_( filter ),
                    upperLimit_( upperLimit ), lowerLimit_( lowerLimit ), maxLevel_( maxLevel ),
                    storage_(storage)
            {}

            void operator()( std::vector< std::pair< const walberla::Block *, uint_t > > & minTargetLevels,
                             std::vector< const walberla::Block * > &,
                             const walberla::BlockForest & )
            {

                for(auto & minTargetLevel : minTargetLevels) {

                    const walberla::Block * const block = minTargetLevel.first;

                    const uint_t currentLevelOfBlock = block->getLevel();

                    const auto * uField = block->template getData< VectorField_T >( fieldID_ );

                    const real_t dx = storage_->dx(currentLevelOfBlock);
                    const real_t dy = storage_->dy(currentLevelOfBlock);
                    const real_t dz = storage_->dz(currentLevelOfBlock);

                    if( uField == nullptr ) {
                        minTargetLevel.second = uint_t(0);
                        continue;
                    }

                    bool refine{ false };
                    bool coarsen{ true };

                    filter_( *block );

                    const auto ci = uField->xyzSizeWithGhostLayer();

                    for( const auto & cell : ci ) {

                        const auto x = cell.x(); const auto y = cell.y(); const auto z = cell.z();

                        real_t qCriterion{};
                        const auto one = cell_idx_t(1);

                        if(filter_(x,y,z) && filter_(x+one,y,z) && filter_(x-one,y,z) && filter_(x,y+one,z)
                           && filter_(x,y-one,z) && filter_(x,y,z+one) && filter_(x,y,z-one)) {
                            const walberla::Vector3<real_t> xa = getVectorFromField(uField,x+one,y,z);
                            const walberla::Vector3<real_t> xb = getVectorFromField(uField,x-one,y,z);
                            const walberla::Vector3<real_t> ya = getVectorFromField(uField,x,y+one,z);
                            const walberla::Vector3<real_t> yb = getVectorFromField(uField,x,y-one,z);
                            const walberla::Vector3<real_t> za = getVectorFromField(uField,x,y,z+one);
                            const walberla::Vector3<real_t> zb = getVectorFromField(uField,x,y,z-one);

                            qCriterion = walberla::lbm::QCriterion::calculate(xa, xb, ya, yb, za, zb, dx, dy, dz);
                        }

                        if( qCriterion > lowerLimit_ ) {
                            coarsen = false;
                            if( qCriterion > upperLimit_ ) {
                                refine = true;
                            }
                        }

                    }

                    if( refine && currentLevelOfBlock < maxLevel_ ) {
                        WALBERLA_ASSERT( !coarsen );
                        minTargetLevel.second = currentLevelOfBlock + uint_t(1);
                    }

                    if( coarsen && currentLevelOfBlock > uint_t(0) ) {
                        WALBERLA_ASSERT( !refine );
                        minTargetLevel.second = currentLevelOfBlock - uint_t(1);
                    }

                }

            } // operator()

        private:

            walberla::Vector3<real_t> getVectorFromField( const VectorField_T * const field,
                                                          const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) {

                walberla::Vector3<real_t> v{};

                if constexpr( VectorField_T::F_SIZE == 1 ) {
                    v = field->get(x,y,z);
                } else if constexpr( VectorField_T::F_SIZE == 3 ) {
                    v[0] = field->get(x,y,z,0);
                    v[1] = field->get(x,y,z,1);
                    v[2] = field->get(x,y,z,2);
                } else {
                    WALBERLA_ABORT("Invalid vector field type.")
                }

                return v;

            }

            const BlockDataID fieldID_;

            Filter_T filter_;

            real_t upperLimit_;
            real_t lowerLimit_;

            uint_t maxLevel_;

            walberla::StructuredBlockStorage * storage_;

        }; // class QCriterionBasedLevelDetermination

    } // namespace refinement

} // namespace turbine_core

#endif // TURBINECORE_QCRITERIONBASEDLEVELDETERMINATION_H 
