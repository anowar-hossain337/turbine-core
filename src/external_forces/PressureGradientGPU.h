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

#ifndef TURBINECORE_PRESSUREGRADIENTGPU_H
#define TURBINECORE_PRESSUREGRADIENTGPU_H

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

#include "external_forces/PressureGradient.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace external_forces {

        namespace internal {
            template<typename Type_T>
            GLOBAL_PREFIX void addPressureGradient( walberla::gpu::FieldAccessor<Type_T> faForce,
                                                    const Type_T pressureGradient, const Type_T windDirection ) {
                faForce.set(blockIdx, threadIdx);
                if(faForce.isValidPosition()) {
                    faForce.get(0) = pressureGradient * cos(windDirection);
                    faForce.get(1) = pressureGradient * sin(windDirection);
                    faForce.get(2) = Type_T(0);
                }
            }
        }

        template< typename GPUField_T >
        struct PressureGradientGPU : public PressureGradientBase<PressureGradientGPU<GPUField_T>, GPUField_T> {

            using Value_T = typename GPUField_T::value_type;
            friend PressureGradientBase<PressureGradientGPU<GPUField_T>, GPUField_T>;

            PressureGradientGPU( const BlockDataID & forceFieldId, const Value_T pressureGradient, const Value_T windDirection,
                                 const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                                 const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
                    : PressureGradientBase<PressureGradientGPU<GPUField_T>, GPUField_T>(
                          forceFieldId, pressureGradient, windDirection, requiredBlockSelectors, incompatibleBlockSelectors
                      )
            {}

        private:

            using FieldIdx_T = walberla::gpu::FieldIndexing<Value_T>;

            virtual void addPressureGradient( GPUField_T * const forceField ) {

                const auto ghostLayer = forceField->nrOfGhostLayers();

                auto pressureKernel = walberla::gpu::make_kernel( &internal::addPressureGradient<Value_T> );
                pressureKernel.addFieldIndexingParam( FieldIdx_T::withGhostLayerXYZ( *forceField, ghostLayer ) );
                pressureKernel.template addParam<Value_T>(this->pressureGradient_);
                pressureKernel.template addParam<Value_T>(this->windDirection_);

                pressureKernel();

            }

        };

    } // namespace external_forces

} // namespace turbine_core

#endif // TURBINECORE_PRESSUREGRADIENTGPU_H