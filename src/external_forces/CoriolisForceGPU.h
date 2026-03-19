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
//! \file CoriolisForceGPU.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//! \author Frédéric BLONDEL <frederic.blondel@ifpen.fr>
//
//======================================================================================================================

#pragma once

#ifndef TURBINECORE_CORIOLISFORCEGPU_H
#define TURBINECORE_CORIOLISFORCEGPU_H

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

#include "external_forces/CoriolisForce.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace external_forces {

        namespace internal {
            template<typename Type_T>
            GLOBAL_PREFIX void addCoriolisForce( walberla::gpu::FieldAccessor<Type_T> faVelocity,
                                                 walberla::gpu::FieldAccessor<Type_T> faForce,
                                                 const Type_T coriolisFrequency, Type_T * const geostrophicWind ) {
                faVelocity.set(blockIdx, threadIdx);
                faForce.set(blockIdx, threadIdx);
                if(faVelocity.isValidPosition() && faForce.isValidPosition()) {
                    faForce.get(0) = - coriolisFrequency * ( geostrophicWind[1] - faVelocity.get(1) );
                    faForce.get(1) = + coriolisFrequency * ( geostrophicWind[0] - faVelocity.get(0) );
                    faForce.get(2) = Type_T(0);
                }
            }
        }

        template< typename GPUField_T >
        struct CoriolisForceGPU : public CoriolisForceBase<CoriolisForceGPU<GPUField_T>, GPUField_T> {

            using Value_T = typename GPUField_T::value_type;
            friend CoriolisForceBase<CoriolisForceGPU<GPUField_T>, GPUField_T>;

            CoriolisForceGPU( const BlockDataID & velocityFieldId, const BlockDataID & forceFieldId,
                              const Value_T coriolisFrequency, const walberla::Vector3<Value_T> & geostrophicWind,
                              const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                              const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
                    : CoriolisForceBase<CoriolisForceGPU<GPUField_T>, GPUField_T>(
                      velocityFieldId, forceFieldId, coriolisFrequency, geostrophicWind,
                      requiredBlockSelectors, incompatibleBlockSelectors)
            {
                WALBERLA_GPU_CHECK(gpuMalloc((void**)&d_geostrophicWind_, 3 * sizeof(Value_T)))
                WALBERLA_GPU_CHECK(gpuMemcpy(d_geostrophicWind_, this->geostrophicWind_.data(), 3 * sizeof(Value_T), gpuMemcpyHostToDevice))
            }

            ~CoriolisForceGPU() { WALBERLA_GPU_CHECK(gpuFree(d_geostrophicWind_)) }

        private:

            using FieldIdx_T = walberla::gpu::FieldIndexing<Value_T>;

            Value_T * const d_geostrophicWind_{};

            virtual void addCoriolisForce( GPUField_T * const velocityField, GPUField_T * const forceField ) {

                const auto ghostLayer = velocityField->nrOfGhostLayers();

                auto coriolisKernel = walberla::gpu::make_kernel( &internal::addCoriolisForce<Value_T> );
                coriolisKernel.addFieldIndexingParam( FieldIdx_T::withGhostLayerXYZ( *velocityField, ghostLayer ) );
                coriolisKernel.addFieldIndexingParam( FieldIdx_T::withGhostLayerXYZ( *forceField, ghostLayer ) );
                coriolisKernel.template addParam<Value_T>(this->coriolisFrequency_);
                coriolisKernel.template addParam<const Value_T*>(d_geostrophicWind_);

                coriolisKernel();

            }

        };

    } // namespace temperature

} // namespace turbine_core

#endif // TURBINECORE_CORIOLISFORCEGPU_H