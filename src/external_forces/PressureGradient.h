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
//! \file PressureGradient.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef TURBINECORE_PRESSUREGRADIENT_H
#define TURBINECORE_PRESSUREGRADIENT_H

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
        struct PressureGradientBase : public Driver {

            using Value_T = typename VectorField_T::value_type;

            PressureGradientBase( const BlockDataID & forceFieldId, const Value_T pressureGradient, const Value_T windDirection,
                                  const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                                  const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
                               : forceFieldId_(forceFieldId), pressureGradient_(pressureGradient), windDirection_(windDirection),
                                 requiredBlockSelectors_(requiredBlockSelectors),
                                 incompatibleBlockSelectors_(incompatibleBlockSelectors)
            {
                WALBERLA_CHECK(!std::isnan(pressureGradient), "You must provide all necessary parameters for the pressure gradient in the input file. \n"
                                                              "These are: initialVelocity and roughnessLengthRatio.")
            }

            void driveFlow(walberla::IBlock * block, const uint_t level, const uint_t executionCount) override {
                this->operator()(block, level, executionCount);
            }

            void operator()(walberla::IBlock * block, const uint_t level = 0, const uint_t executionCount = 0)
            {
                if( !walberla::selectable::isSetSelected( block->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_ ) ) {
                    return;
                }

                auto * forceField = block->getData<VectorField_T>(forceFieldId_);
                WALBERLA_ASSERT_NOT_NULLPTR(forceField, "Invalid force field.")
                WALBERLA_ASSERT(forceField->fSize() == 3, "You must provide a three-dimensional force field.")

                static_cast<T*>(this)->addPressureGradient(forceField);

            }

        protected:

            const Value_T pressureGradient_;
            const Value_T windDirection_;

        private:

            virtual void addPressureGradient( VectorField_T * const forceField ) = 0;

            const BlockDataID forceFieldId_;

            const walberla::Set<walberla::SUID> requiredBlockSelectors_;
            const walberla::Set<walberla::SUID> incompatibleBlockSelectors_;
        };

        template< typename VectorField_T >
        struct PressureGradient : public PressureGradientBase<PressureGradient<VectorField_T>, VectorField_T> {

            using Value_T = typename VectorField_T::value_type;
            friend PressureGradientBase<PressureGradient<VectorField_T>, VectorField_T>;

            PressureGradient( const BlockDataID & forceFieldId, const Value_T pressureGradient, const Value_T windDirection,
                              const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                              const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
                           : PressureGradientBase<PressureGradient<VectorField_T>, VectorField_T>(
                                 forceFieldId, pressureGradient, windDirection, requiredBlockSelectors, incompatibleBlockSelectors
                            )
            {}

        private:

            virtual void addPressureGradient( VectorField_T * const forceField ) {
                for( auto cellIt = forceField->beginWithGhostLayer(); cellIt != forceField->end(); ++cellIt) {
                    // Assumes the flow is driven in the x-direction
                    forceField->get(cellIt.cell(), 0) = this->pressureGradient_ * cos(this->windDirection_);
                    forceField->get(cellIt.cell(), 1) = this->pressureGradient_ * sin(this->windDirection_);
                    forceField->get(cellIt.cell(), 2) = Value_T(0);
                }
            }

        };

    } // namespace external_forces

} // namespace turbine_core

#endif // TURBINECORE_PRESSUREGRADIENT_H