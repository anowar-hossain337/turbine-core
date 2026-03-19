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
//! \file PressureGradientGPU.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef TURBINECORE_DYNAMICPRESSUREGRADIENTGPU_H
#define TURBINECORE_DYNAMICPRESSUREGRADIENTGPU_H

#include <functional>

#include <core/Set.h>
#include <core/selectable/IsSetSelected.h>

#include <core/uid/SUID.h>

#include <domain_decomposition/IBlock.h>
#include <blockforest/BlockForest.h>
#include <blockforest/StructuredBlockForest.h>
#include <gpu/Kernel.h>
#include <gpu/FieldAccessor.h>
#include <gpu/FieldIndexing.h>

#include "external_forces/DynamicPressureGradient.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace external_forces {

        namespace internal {
            template<typename Type_T>
            GLOBAL_PREFIX void addDynamicPressureGradient( walberla::gpu::FieldAccessor<Type_T> faVelocity,
                                                    walberla::gpu::FieldAccessor<Type_T> faForce,
                                                    const Type_T pressureGradient, const Type_T windDirection, const Type_T windDirectionChangeRate ) {
                faVelocity.set(blockIdx, threadIdx);
                faForce.set(blockIdx, threadIdx);
                if(faVelocity.isValidPosition() && faForce.isValidPosition()) {
                    faForce.get(0) = pressureGradient * cos(windDirection) - real_t(0.5) * windDirectionChangeRate * faVelocity.get(1);
                    faForce.get(1) = pressureGradient * sin(windDirection) + real_t(0.5) * windDirectionChangeRate * faVelocity.get(0);
                    faForce.get(2) = Type_T(0);
                }
            }
        }

        template< typename GPUField_T >
        struct DynamicPressureGradientGPU : public DynamicPressureGradientBase<DynamicPressureGradientGPU<GPUField_T>, GPUField_T> {

            using Value_T = typename GPUField_T::value_type;
            friend DynamicPressureGradientBase<DynamicPressureGradientGPU<GPUField_T>, GPUField_T>;

            DynamicPressureGradientGPU( const BlockDataID & velocityFieldId, const BlockDataID & forceFieldId, const Value_T pressureGradient, const Value_T windDirection,
                                walberla::timeloop::ITimeloop * timeloop, const data_interpolation::DataInterpolator<real_t, real_t, false> windDirChangeRateInterpolator,
                                 const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                                 const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
                    : DynamicPressureGradientBase<DynamicPressureGradientGPU<GPUField_T>, GPUField_T>(
                          velocityFieldId, forceFieldId, pressureGradient, windDirection, timeloop, windDirChangeRateInterpolator, requiredBlockSelectors, incompatibleBlockSelectors
                      ) {}

        private:

            using FieldIdx_T = walberla::gpu::FieldIndexing<Value_T>;

            virtual void addDynamicPressureGradient( GPUField_T * const velocityField, GPUField_T * const forceField ) {

                const auto ghostLayer = velocityField->nrOfGhostLayers();

                auto pressureKernel = walberla::gpu::make_kernel( &internal::addDynamicPressureGradient<Value_T> );
                pressureKernel.addFieldIndexingParam( FieldIdx_T::withGhostLayerXYZ( *velocityField, ghostLayer ) );
                pressureKernel.addFieldIndexingParam( FieldIdx_T::withGhostLayerXYZ( *forceField, ghostLayer ) );
                pressureKernel.template addParam<Value_T>(this->pressureGradient_);
                pressureKernel.template addParam<Value_T>(this->windDirection_);
                pressureKernel.template addParam<Value_T>(this->windDirChangeRate_);

                pressureKernel();

            }

        };

    } // namespace external_forces

} // namespace turbine_core

#endif // TURBINECORE_DYNAMICPRESSUREGRADIENTGPU_H