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
//! \file FieldCopy.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef TURBINECORE_FIELDCOPY_H
#define TURBINECORE_FIELDCOPY_H

#include <domain_decomposition/StructuredBlockStorage.h>
#include <field/Layout.h>

#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace field {

        template<typename Field_T>
        void fieldCopy(Field_T * dst, const Field_T * src ) {

            if (dst->layout() != src->layout()) { WALBERLA_ABORT("Cannot copy fields with different layout") }

            bool canCopy =
                    (src->layout() == walberla::field::fzyx && dst->fAllocSize() == src->fAllocSize() && dst->zAllocSize() == src->zAllocSize() &&
                     dst->yAllocSize() == src->yAllocSize() && dst->xSize() == src->xSize()) ||
                    (src->layout() == walberla::field::zyxf && dst->zAllocSize() == src->zAllocSize() && dst->yAllocSize() == src->yAllocSize() &&
                     dst->xAllocSize() == src->xAllocSize() && dst->fSize() == src->fSize());

            if (!canCopy) { WALBERLA_ABORT("Field have to have the same size ") }

            std::copy( src->data(), src->data() + src->allocSize(), dst->data() );

        }

        template<typename Field_T>
        void fieldCopy( const std::shared_ptr< walberla::StructuredBlockStorage > & blocks,
                        BlockDataID dstID, const BlockDataID srcID )
        {
            for ( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
            {
                auto dst = blockIt->getData<Field_T>( dstID );
                const auto src = blockIt->getData<Field_T>( srcID );
                fieldCopy( dst, src );
            }
        }

    } // namespace gpu

} // namespace turbine_core

#endif // TURBINECORE_FIELDCOPY_H