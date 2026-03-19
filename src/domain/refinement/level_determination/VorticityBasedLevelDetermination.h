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
//! \file VorticityBasedLevelDetermination.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_VORTICITYBASEDLEVELDETERMINATION_H
#define TURBINECORE_VORTICITYBASEDLEVELDETERMINATION_H

#pragma once

#include <utility>
#include <vector>

#include <blockforest/Block.h>
#include <core/math/Vector3.h>

#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace refinement {

        /*
         * Refinement check based on vorticity
         * If vorticity magnitude is below lowerLimit in all cells of a block, that block could be coarsened.
         * If the vorticity is above the upperLimit for at least one cell, that block gets marked for refinement.
         * Else, the block remains on the current level.
         */
        template< typename Stencil_T, typename VectorField_T, typename Filter_T >
        class VorticityBasedLevelDetermination
        {
        public:

            VorticityBasedLevelDetermination( const BlockDataID & fieldID, const Filter_T & filter,
                                                   const real_t upperLimit, const real_t lowerLimit, const uint_t maxLevel )
            : fieldID_( fieldID ), filter_( filter ),
              upperLimitSqr_( upperLimit*upperLimit ), lowerLimitSqr_( lowerLimit*lowerLimit ), maxLevel_( maxLevel )
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

                    walberla::Vector3<real_t> uCurl{};

                    bool refine{ false };
                    bool coarsen{ true };

                    filter_( *block );

                    const auto ci = uField->xyzSizeWithGhostLayer();

                    for( const auto & cell : ci ) {

                        const auto x = cell.x(); const auto y = cell.y(); const auto z = cell.z();

                        std::vector< walberla::Vector3<real_t> > uValues( Stencil_T::Size, walberla::Vector3<real_t>{} );
                        getVelocities(uValues,x,y,z,uField);

                        // obtain the vector curl(u) with the help of the gradient formula from
                        // See: Ramadugu et al - Lattice differential operators for computational physics (2013)
                        // higher order terms are neglected
                        // with T = c_s**2
                        const real_t inv_c_s_sqr{3};
                        uCurl.reset();
                        calculateCurl(uCurl, uValues);

                        uCurl *= inv_c_s_sqr;

                        // squared 2 norm of curl
                        real_t normSq = uCurl.sqrLength();

                        if( normSq > lowerLimitSqr_ ) {
                            coarsen = false;
                            if( normSq > upperLimitSqr_ ) {
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

            real_t upperLimitSqr_;
            real_t lowerLimitSqr_;

            uint_t maxLevel_;

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

            void calculateCurl( walberla::Vector3<real_t> & curl, const std::vector<walberla::Vector3<real_t>> & velocities ) {

                for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir) {

                    const walberla::Vector3<real_t> latticeVelocity{ real_t(dir.cx()), real_t(dir.cy()), real_t(dir.cz()) };
                    const walberla::Vector3<real_t> & velocity = velocities[*dir];

                    const real_t w_i = LatticeModel_T::w[dir.toIdx()];

                    const walberla::Vector3<real_t> curl_i = latticeVelocity % velocity;

                    curl += w_i * curl_i;

                }

            }

        }; // class GradientRefinement

    } // namespace refinement

} // namespace turbine_core

#endif // TURBINECORE_VORTICITYBASEDLEVELDETERMINATION_H 
