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
//! \file VectorGradientBasedLevelDetermination.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_VECTORGRADIENTBASEDLEVELDETERMINATION_H
#define TURBINECORE_VECTORGRADIENTBASEDLEVELDETERMINATION_H

#pragma once

#include <utility>
#include <vector>

#include <blockforest/Block.h>
#include <core/math/Vector3.h>
#include <core/math/Matrix3.h>

#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace refinement {

        /*
         * Refinement check based on gradient magnitude
         * If gradient magnitude is below lowerLimit in all cells of a block, that block could be coarsened.
         * If the gradient value is above the upperLimit for at least one cell, that block gets marked for refinement.
         * Else, the block remains on the current level.
         *
         * DISCLAIMER : heavily based on Christoph R's version.
         */
        template< typename Stencil_T, typename VectorField_T, typename Filter_T >
        class VectorGradientBasedLevelDetermination
        {
        public:

            VectorGradientBasedLevelDetermination( const BlockDataID & fieldID, const Filter_T & filter, walberla::StructuredBlockStorage * const storage,
                                                   const real_t upperLimit, const real_t lowerLimit, const uint_t maxLevel ) :
                    fieldID_( fieldID ), filter_( filter ),
                    upperLimit_( upperLimit ), lowerLimit_( lowerLimit ), maxLevel_( maxLevel ), storage_(storage)
            {}

            void operator()( std::vector< std::pair< const walberla::Block *, uint_t > > & minTargetLevels,
                             std::vector< const walberla::Block * > &,
                             const walberla::BlockForest & )         {

                for(auto & minTargetLevel : minTargetLevels) {

                    const walberla::Block * const block = minTargetLevel.first;

                    const uint_t currentLevelOfBlock = block->getLevel();

                    const auto * uField = block->template getData< VectorField_T >( fieldID_ );

                    if( uField == nullptr ) {
                        minTargetLevel.second = uint_t(0);
                        continue;
                    }

                    walberla::Matrix3<real_t> uGradient( real_t(0) );

                    bool refine{ false };
                    bool coarsen{ true };

                    filter_( *block );

                    const auto ci = uField->xyzSizeWithGhostLayer();

                    for( const auto & cell : ci ) {

                        const cell_idx_t x = cell.x(); const cell_idx_t y = cell.y(); const cell_idx_t z = cell.z();

                        std::vector< walberla::Vector3<real_t> > uValues( Stencil_T::Size, walberla::Vector3<real_t>{} );
                        getVelocities(uValues,x,y,z,uField);

                        // obtain the matrix grad(u) with the help of the gradient formula from
                        // See: Ramadugu et al - Lattice differential operators for computational physics (2013)
                        // with T = c_s**2
                        const real_t inv_c_s_sqr{3};
                        uGradient = real_t(0);
                        calculateGradient(uGradient, uValues);

                        uGradient *= inv_c_s_sqr;

                        real_t norm{};
                        //compute maximums norm of 3x3 matrix
                        for( uint_t i = 0; i < 9; ++i ) {
                            norm = std::max(norm, std::fabs(uGradient[i]));
                        }

                        if( norm > lowerLimit_ ) {
                            coarsen = false;
                            if( norm > upperLimit_ ) {
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

            static constexpr uint_t fSize_ = VectorField_T::F_SIZE;

            const BlockDataID fieldID_;

            Filter_T filter_;

            real_t upperLimit_;
            real_t lowerLimit_;

            uint_t maxLevel_;

            walberla::StructuredBlockStorage * const storage_{};

            void getVelocities( std::vector<walberla::Vector3<real_t>> & velocities, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const VectorField_T * uField ) {

                walberla::Vector3<real_t> velocityInCellCenter{};

                if constexpr (fSize_ == 1) {
                    velocityInCellCenter = uField->get(x,y,z);
                } else if constexpr (fSize_ == 3) {
                    velocityInCellCenter[0] = uField->get( x,y,z,0 );
                    velocityInCellCenter[1] = uField->get( x,y,z,1 );
                    velocityInCellCenter[2] = uField->get( x,y,z,2 );
                } else {
                    WALBERLA_ABORT("Invalid fSize for Vector Field.")
                }

                for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir) {
                    // check if boundary treatment is necessary
                    if( filter_( x+dir.cx(),y+dir.cy(),z+dir.cz() ) ) {
                        // copy from center cell
                        velocities[*dir] = velocityInCellCenter;
                    } else {
                        if constexpr (fSize_ == 1) {
                            velocities[*dir] = uField->get( x+dir.cx(),y+dir.cy(),z+dir.cz() );
                        } else if constexpr (fSize_ == 3) {
                            velocities[*dir][0] = uField->get( x+dir.cx(),y+dir.cy(),z+dir.cz(),0 );
                            velocities[*dir][1] = uField->get( x+dir.cx(),y+dir.cy(),z+dir.cz(),1 );
                            velocities[*dir][2] = uField->get( x+dir.cx(),y+dir.cy(),z+dir.cz(),2 );
                        } else {
                            WALBERLA_ABORT("Invalid fSize for Vector Field.")
                        }
                    }
                }
            }

            void calculateGradient( walberla::Matrix3<real_t> & gradient, const std::vector<walberla::Vector3<real_t>> & velocities ) {

                for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir) {

                    const real_t cx = real_t(dir.cx());
                    const real_t cy = real_t(dir.cy());
                    const real_t cz = real_t(dir.cz());

                    // grad(ux)
                    const real_t ux = velocities[*dir][0];
                    gradient[0] += LatticeModel_T::w[dir.toIdx()] * cx * ux;
                    gradient[3] += LatticeModel_T::w[dir.toIdx()] * cy * ux;
                    gradient[6] += LatticeModel_T::w[dir.toIdx()] * cz * ux;

                    // grad(uy)
                    const real_t uy = velocities[*dir][1];
                    gradient[1] += LatticeModel_T::w[dir.toIdx()] * cx * uy;
                    gradient[4] += LatticeModel_T::w[dir.toIdx()] * cy * uy;
                    gradient[7] += LatticeModel_T::w[dir.toIdx()] * cz * uy;

                    // grad(uz)
                    const real_t uz = velocities[*dir][2];
                    gradient[2] += LatticeModel_T::w[dir.toIdx()] * cx * uz;
                    gradient[5] += LatticeModel_T::w[dir.toIdx()] * cy * uz;
                    gradient[8] += LatticeModel_T::w[dir.toIdx()] * cz * uz;

                }

            }

        }; // class VectorGradientBasedLevelDetermination

    } // namespace refinement

} // namespace turbine_core

#endif // TURBINECORE_VECTORGRADIENTBASEDLEVELDETERMINATION_H 
