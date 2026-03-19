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
//! \file FieldCopyGPU.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef TURBINECORE_FIELDCOPY_GPU_H
#define TURBINECORE_FIELDCOPY_GPU_H

#include <domain_decomposition/StructuredBlockStorage.h>
#include <gpu/GPUField.h>

#include "walberla_helper/field/FieldCopy.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace gpu {

        template<typename T>
        void fieldCopy(walberla::gpu::GPUField<T> * dst, const walberla::gpu::GPUField<T> * src ) {

            gpuMemcpy3DParms p;
            memset(&p, 0, sizeof(p));

            if (dst->layout() != src->layout()) { WALBERLA_ABORT("Cannot copy fields with different layout") }

            bool canCopy =
                    (src->layout() == walberla::field::fzyx && dst->fAllocSize() == src->fAllocSize() && dst->zAllocSize() == src->zAllocSize() &&
                     dst->yAllocSize() == src->yAllocSize() && dst->xSize() == src->xSize()) ||
                    (src->layout() == walberla::field::zyxf && dst->zAllocSize() == src->zAllocSize() && dst->yAllocSize() == src->yAllocSize() &&
                     dst->xAllocSize() == src->xAllocSize() && dst->fSize() == src->fSize());

            if (!canCopy) { WALBERLA_ABORT("Field have to have the same size ") }

            p.srcPtr = src->pitchedPtr();
            p.dstPtr = dst->pitchedPtr();
            p.kind   = gpuMemcpyDeviceToDevice;
            WALBERLA_GPU_CHECK(gpuMemcpy3D(&p))

        }

    } // namespace gpu

} // namespace turbine_core

#endif // TURBINECORE_FIELDCOPY_GPU_H