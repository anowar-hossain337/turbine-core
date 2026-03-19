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
//! \file DynamicPressureGradient.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef TURBINECORE_DYNAMICPRESSUREGRADIENT_H
#define TURBINECORE_DYNAMICPRESSUREGRADIENT_H

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
        struct DynamicPressureGradientBase : public Driver {

            using Value_T = typename VectorField_T::value_type;

            DynamicPressureGradientBase( const BlockDataID & velocityFieldId, const BlockDataID & forceFieldId, const Value_T pressureGradient, const Value_T windDirection,
                                  walberla::timeloop::ITimeloop * timeloop, const data_interpolation::DataInterpolator<real_t, real_t, false> windDirChangeRateInterpolator,
                                  const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                                  const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
                               : velocityFieldId_(velocityFieldId), forceFieldId_(forceFieldId), pressureGradient_(pressureGradient), windDirection_(windDirection),
                                 timeloop_(timeloop), windDirChangeRateInterpolator_(windDirChangeRateInterpolator),
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

                auto * velocityField = block->getData<VectorField_T>(velocityFieldId_);
                WALBERLA_ASSERT_NOT_NULLPTR(velocityField, "Invalid velocity field.")
                WALBERLA_ASSERT(velocityField->fSize() == 3, "You must provide a three-dimensional velocity field.")

                auto * forceField = block->getData<VectorField_T>(forceFieldId_);
                WALBERLA_ASSERT_NOT_NULLPTR(forceField, "Invalid force field.")
                WALBERLA_ASSERT(forceField->fSize() == 3, "You must provide a three-dimensional force field.")

                // update the wind direction
                windDirChangeRate_ = windDirChangeRateInterpolator_(timeloop_->getCurrentTimeStep());
                windDirection_ += windDirChangeRate_;

                static_cast<T*>(this)->addDynamicPressureGradient(velocityField, forceField);

            }

        protected:

            const Value_T pressureGradient_;
            Value_T windDirection_;
            Value_T windDirChangeRate_;

            walberla::timeloop::ITimeloop * timeloop_{nullptr};
            const data_interpolation::DataInterpolator<real_t, real_t, false> windDirChangeRateInterpolator_;

        private:

            virtual void addDynamicPressureGradient( VectorField_T * const velocityField, VectorField_T * const forceField ) = 0;

            const BlockDataID velocityFieldId_;
            const BlockDataID forceFieldId_;

            const walberla::Set<walberla::SUID> requiredBlockSelectors_;
            const walberla::Set<walberla::SUID> incompatibleBlockSelectors_;
        };

        template< typename VectorField_T >
        struct DynamicPressureGradient : public DynamicPressureGradientBase<DynamicPressureGradient<VectorField_T>, VectorField_T> {

            using Value_T = typename VectorField_T::value_type;
            friend DynamicPressureGradientBase<DynamicPressureGradient<VectorField_T>, VectorField_T>;

            DynamicPressureGradient( const BlockDataID & velocityFieldId, const BlockDataID & forceFieldId, const Value_T pressureGradient, const Value_T windDirection,
                                     walberla::timeloop::ITimeloop * timeloop, const data_interpolation::DataInterpolator<real_t, real_t, false> windDirChangeRateInterpolator,
                                     const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                                     const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
                           : DynamicPressureGradientBase<DynamicPressureGradient<VectorField_T>, VectorField_T>(
                                 velocityFieldId, forceFieldId, pressureGradient, windDirection, timeloop, windDirChangeRateInterpolator, requiredBlockSelectors, incompatibleBlockSelectors
                            )
            {}

        private:

            virtual void addDynamicPressureGradient( VectorField_T * const velocityField,  VectorField_T * const forceField ) {
                for( auto cellIt = forceField->beginWithGhostLayer(); cellIt != forceField->end(); ++cellIt) {
                    // Assumes the flow is driven in the x-direction
                    forceField->get(cellIt.cell(), 0) = this->pressureGradient_ * cos(this->windDirection_) - real_t(0.5) * this->windDirChangeRate_ * velocityField->get(cellIt.cell(), 1);
                    forceField->get(cellIt.cell(), 1) = this->pressureGradient_ * sin(this->windDirection_) + real_t(0.5) * this->windDirChangeRate_ * velocityField->get(cellIt.cell(), 0);
                    forceField->get(cellIt.cell(), 2) = Value_T(0);
                }
            }

        };

    } // namespace external_forces

} // namespace turbine_core

#endif // TURBINECORE_DYNAMICPRESSUREGRADIENT_H