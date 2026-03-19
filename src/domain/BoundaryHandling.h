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
//! \file BoundaryHandling.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_BOUNDARYHANDLING_H
#define TURBINECORE_BOUNDARYHANDLING_H

#pragma once

#include <blockforest/BlockDataHandling.h>

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/ProjectDefines.h"
#include "EnvironmentSetup.h"

#include "wind_turbine_core/math/Vector3.h"

namespace turbine_core {

    namespace boundary {
        namespace internal {
            template< typename Stencil_T, typename FlagField_T, typename... Boundaries >
            using BoundaryHandling_T = walberla::BoundaryHandling<FlagField_T, Stencil_T, Boundaries...>;

            template< typename Stencil_T, typename FlagField_T, typename Boundary_T>
            using DummyBoundaryHandling_T = walberla::BoundaryHandling<FlagField_T, Stencil_T, Boundary_T>;
        }

        auto VelocityCallback = [](const walberla::Cell& cellCenter,
                                   const std::shared_ptr<walberla::StructuredBlockForest >& blockForest,
                                   const walberla::IBlock& block, const real_t roughnessLength, const real_t kappa, const real_t frictionVelocity) {

            walberla::CellInterval domain = blockForest->getDomainCellBB();

            walberla::Cell globalCell;
            blockForest->transformBlockLocalToGlobalCell(globalCell, block, cellCenter);

            walberla::Vector3<real_t> localCellCenter;
            blockForest->getCellCenter(localCellCenter, globalCell, blockForest->getLevel(block));

            const real_t height = walberla::math::max(localCellCenter[2], real_t(0.05));
            auto velocity = frictionVelocity / kappa * std::log( height / roughnessLength );

            walberla::math::Vector3< real_t > result(velocity, 0.0, 0.0);

            return result;
        };

        auto velocityInit(const real_t roughnessLength, const real_t kappa, const real_t frictionVelocity) {
            std::function< walberla::math::Vector3< real_t >(const walberla::Cell&, const std::shared_ptr<walberla::StructuredBlockForest >&, walberla::IBlock&)>
                    velocityInitialisation = std::bind(turbine_core::boundary::VelocityCallback, std::placeholders::_1,
                                                       std::placeholders::_2, std::placeholders::_3, roughnessLength, kappa, frictionVelocity);
            return velocityInitialisation;
        }

        template< typename Stencil_T, typename PdfField_T, typename FlagField_T, typename DummyNoSlip_T>
        class DummyBoundaryHandling
                : public walberla::blockforest::AlwaysInitializeBlockDataHandling<internal::DummyBoundaryHandling_T<Stencil_T, FlagField_T, DummyNoSlip_T>> {

        public:

            DummyBoundaryHandling(const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID)
                    : flagFieldID_(flagFieldID), pdfFieldID_(pdfFieldID) {}

            using BoundaryHandling = internal::DummyBoundaryHandling_T<Stencil_T, FlagField_T, DummyNoSlip_T>;

            BoundaryHandling * initialize(walberla::IBlock * const block) {
                auto * flagField = block->getData<FlagField_T>(flagFieldID_);
                auto * pdfField = block->getData<PdfField_T>(pdfFieldID_);
                const auto fluidFlag = flagField->getOrRegisterFlag("DummyFlag");

                auto * handling = new BoundaryHandling(
                        "Dummy Boundary Handling", flagField, fluidFlag,
                        DummyNoSlip_T("DummyNoSlip", walberla::FlagUID("DummyWall"), pdfField));

                return handling;
            }

        private:
            const BlockDataID flagFieldID_;
            const BlockDataID pdfFieldID_;

        }; // class TurbineBoundaryHandling
    }
}

#endif