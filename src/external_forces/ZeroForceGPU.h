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
//! \file ZeroForceGPU.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef TURBINECORE_ZEROFORCEGPU_H
#define TURBINECORE_ZEROFORCEGPU_H

#include <functional>

#include <core/Set.h>
#include <core/selectable/IsSetSelected.h>

#include <core/uid/SUID.h>

#include "external_forces/ZeroForce.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace external_forces {

        template< typename GPUField_T >
        struct ZeroForceGPU : public ZeroForceBase<ZeroForceGPU<GPUField_T>, GPUField_T> {

            using Value_T = typename GPUField_T::value_type;
            friend ZeroForceBase<ZeroForceGPU<GPUField_T>, GPUField_T>;

            ZeroForceGPU( const BlockDataID & forceFieldId,
                          const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                          const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
                    : ZeroForceBase<ZeroForceGPU<GPUField_T>, GPUField_T>(
                          forceFieldId, requiredBlockSelectors, incompatibleBlockSelectors
                      )
            {}

        private:

            virtual void addZeroForce( GPUField_T * const forceField ) {
                WALBERLA_GPU_CHECK(gpuMemset(forceField->data(), Value_T(0), forceField->allocSize() * sizeof(Value_T)))
            }

        };

    } // namespace external_forces

} // namespace turbine_core

#endif // TURBINECORE_ZEROFORCEGPU_H