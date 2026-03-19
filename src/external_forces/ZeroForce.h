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
//! \file ZeroForce.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef TURBINECORE_ZEROFORCE_H
#define TURBINECORE_ZEROFORCE_H

#include <functional>

#include <core/Set.h>
#include <core/selectable/IsSetSelected.h>

#include <core/uid/SUID.h>

#include <domain_decomposition/IBlock.h>

#include "external_forces/DriverSetup.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace external_forces {

        template< typename T, typename VectorField_T >
        struct ZeroForceBase : public Driver {

            using Value_T = typename VectorField_T::value_type;

            ZeroForceBase( const BlockDataID & forceFieldId,
                           const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                           const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
                           : forceFieldId_(forceFieldId),
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

                auto * forceField = block->getData<VectorField_T>(forceFieldId_);
                WALBERLA_ASSERT_NOT_NULLPTR(forceField, "Invalid force field.")

                static_cast<T*>(this)->addZeroForce(forceField);

            }

        private:

            virtual void addZeroForce( VectorField_T * const forceField ) = 0;

            const BlockDataID forceFieldId_;

            const walberla::Set<walberla::SUID> requiredBlockSelectors_;
            const walberla::Set<walberla::SUID> incompatibleBlockSelectors_;
        };

        template< typename VectorField_T >
        struct ZeroForce : public ZeroForceBase<ZeroForce<VectorField_T>, VectorField_T> {

            using Value_T = typename VectorField_T::value_type;
            friend ZeroForceBase<ZeroForce<VectorField_T>, VectorField_T>;

            ZeroForce( const BlockDataID & forceFieldId,
                  const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                  const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
                : ZeroForceBase<ZeroForce<VectorField_T>, VectorField_T>(
                      forceFieldId, requiredBlockSelectors, incompatibleBlockSelectors
                 )
            {}

        private:

            virtual void addZeroForce( VectorField_T * const forceField ) {
                forceField->setWithGhostLayer(Value_T(0));
            }

        };

    } // namespace external_forces

} // namespace turbine_core

#endif // TURBINECORE_ZEROFORCE_H