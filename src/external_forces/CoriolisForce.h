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
//! \file CoriolisForce.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//! \author Frédéric BLONDEL <frederic.blondel@ifpen.fr>
//
//======================================================================================================================

#pragma once

#ifndef TURBINECORE_CORIOLISFORCE_H
#define TURBINECORE_CORIOLISFORCE_H

#include <functional>

#include <core/Set.h>
#include <core/selectable/IsSetSelected.h>

#include <core/uid/SUID.h>

#include <domain_decomposition/IBlock.h>
#include <blockforest/BlockForest.h>
#include <blockforest/StructuredBlockForest.h>

#include "external_forces/DriverSetup.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace external_forces {

        template< typename T, typename VectorField_T >
        struct CoriolisForceBase : public Driver {

            using Value_T = typename VectorField_T::value_type;

            CoriolisForceBase( const BlockDataID & velocityFieldId, const BlockDataID & forceFieldId,
                               const Value_T coriolisFrequency, const walberla::Vector3<Value_T> & geostrophicWind,
                               const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                               const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
                               : velocityFieldId_(velocityFieldId), forceFieldId_(forceFieldId),
                                 coriolisFrequency_(coriolisFrequency), geostrophicWind_(geostrophicWind),
                                 requiredBlockSelectors_(requiredBlockSelectors),
                                 incompatibleBlockSelectors_(incompatibleBlockSelectors)
            {}

            void driveFlow(walberla::IBlock * block, const uint_t level, const uint_t executionCount) override {
                this->operator()(block, level, executionCount);
            }

            void operator()(walberla::IBlock * block, const uint_t level = 0, const uint_t executionCount = 0)
            {
                if( !walberla::selectable::isSetSelected( block->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_ ) ) {
                    return;
                }

                auto * velocityField = block->getData<VectorField_T>(velocityFieldId_);
                WALBERLA_ASSERT_NOT_NULLPTR(velocityField, "Invalid velocity field.")
                WALBERLA_ASSERT(velocityField->fSize() == 3, "You must provide a three-dimensional velocity field.")

                auto * forceField = block->getData<VectorField_T>(forceFieldId_);
                WALBERLA_ASSERT_NOT_NULLPTR(forceField, "Invalid force field.")
                WALBERLA_ASSERT(forceField->fSize() == 3, "You must provide a three-dimensional force field.")

                static_cast<T*>(this)->addCoriolisForce(velocityField, forceField);

            }

        protected:

            const Value_T coriolisFrequency_;
            const walberla::Vector3<Value_T> geostrophicWind_;

        private:

            virtual void addCoriolisForce( VectorField_T * const velocityField, VectorField_T * const forceField ) = 0;

            const BlockDataID velocityFieldId_;
            const BlockDataID forceFieldId_;

            const walberla::Set<walberla::SUID> requiredBlockSelectors_;
            const walberla::Set<walberla::SUID> incompatibleBlockSelectors_;
        };

        template< typename VectorField_T >
        struct CoriolisForce : public CoriolisForceBase<CoriolisForce<VectorField_T>, VectorField_T> {

            using Value_T = typename VectorField_T::value_type;
            friend CoriolisForceBase<CoriolisForce<VectorField_T>, VectorField_T>;

            CoriolisForce( const BlockDataID & velocityFieldId, const BlockDataID & forceFieldId,
                           const Value_T coriolisFrequency, const walberla::Vector3<Value_T> & geostrophicWind,
                           const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                           const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
                           : CoriolisForceBase<CoriolisForce<VectorField_T>, VectorField_T>(
                                 velocityFieldId, forceFieldId, coriolisFrequency, geostrophicWind,
                                 requiredBlockSelectors, incompatibleBlockSelectors
                            )
            {

                WALBERLA_CHECK(!std::isnan(coriolisFrequency), "You must provide all necessary parameters for the Coriolis frequency in the input file. \n")
                WALBERLA_CHECK(!std::isnan(geostrophicWind[0]), "You must provide all necessary parameters for the geostrophic wind in the input file. \n")
            }

        private:

            virtual void addCoriolisForce( VectorField_T * const velocityField, VectorField_T * const forceField ) {
                for( auto cellIt = forceField->beginWithGhostLayer(); cellIt != forceField->end(); ++cellIt) {

                    // Assumes the third component is in the vertical direction
                    forceField->get(cellIt.cell(), 0) = - this->coriolisFrequency_ * ( this->geostrophicWind_[1] - velocityField->get(cellIt.cell(), 1) );
                    forceField->get(cellIt.cell(), 1) = + this->coriolisFrequency_ * ( this->geostrophicWind_[0] - velocityField->get(cellIt.cell(), 0) );
                    forceField->get(cellIt.cell(), 2) = Value_T(0);
                }
            }

        };

    } // namespace field

} // namespace turbine_core

#endif // TURBINECORE_CORIOLISFORCE_H