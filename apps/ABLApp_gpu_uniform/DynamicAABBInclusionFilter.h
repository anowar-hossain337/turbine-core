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
//! \file DynamicAABBInclusionFilter.h
//
//======================================================================================================================

#ifndef TURBINECORE_DYNAMICAABBINCLUSIONFILTER_H
#define TURBINECORE_DYNAMICAABBINCLUSIONFILTER_H

#pragma once

#include <memory>

#include <core/DataTypes.h>
#include <core/cell/CellSet.h>
#include <core/math/AABB.h>
#include <domain_decomposition/IBlock.h>
#include <domain_decomposition/StructuredBlockStorage.h>

#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {
namespace output {

// Reason for edit start: provide an app-local inclusion filter whose AABB bounds can be updated at runtime without modifying the stock waLBerla AABB cell filter.
class DynamicAABBInclusionFilter {
public:
    using Vector_T = walberla::Vector3<real_t>;

    class State {
    public:
        State(const Vector_T & minCorner, const Vector_T & maxCorner)
                : aabb_(minCorner, maxCorner)
        {}

        explicit State(const walberla::AABB & aabb)
                : aabb_(aabb)
        {}

        void setBounds(const Vector_T & minCorner, const Vector_T & maxCorner)
        {
            aabb_ = walberla::AABB(minCorner, maxCorner);
        }

        void setAABB(const walberla::AABB & aabb)
        {
            aabb_ = aabb;
        }

        const walberla::AABB & getAABB() const
        {
            return aabb_;
        }

    private:
        walberla::AABB aabb_;
    };

    explicit DynamicAABBInclusionFilter(const std::shared_ptr<State> & state)
            : state_(state)
    {
        WALBERLA_ASSERT_NOT_NULLPTR(state_);
    }

    void operator()(walberla::CellSet & filteredCells,
                    const walberla::IBlock & block,
                    const walberla::StructuredBlockStorage & storage,
                    const uint_t ghostLayers = uint_t(0)) const
    {
        WALBERLA_ASSERT_NOT_NULLPTR(state_);

        walberla::CellInterval cellBB;
        storage.getCellBBFromAABB(cellBB, state_->getAABB(), storage.getLevel(block));
        storage.transformGlobalToBlockLocalCellInterval(cellBB, block);

        const walberla::cell_idx_t start = walberla::cell_idx_c(-1) * walberla::cell_idx_c(ghostLayers);

        for(walberla::cell_idx_t z = start; z != walberla::cell_idx_c(storage.getNumberOfZCells(block)) + walberla::cell_idx_c(ghostLayers); ++z)
            for(walberla::cell_idx_t y = start; y != walberla::cell_idx_c(storage.getNumberOfYCells(block)) + walberla::cell_idx_c(ghostLayers); ++y)
                for(walberla::cell_idx_t x = start; x != walberla::cell_idx_c(storage.getNumberOfXCells(block)) + walberla::cell_idx_c(ghostLayers); ++x)
                    if(cellBB.contains(x, y, z))
                        filteredCells.insert(x, y, z);
    }

private:
    std::shared_ptr<State> state_;
};
// Reason for edit end: provide an app-local inclusion filter whose AABB bounds can be updated at runtime without modifying the stock waLBerla AABB cell filter.

} // namespace output
} // namespace turbine_core

#endif // TURBINECORE_DYNAMICAABBINCLUSIONFILTER_H
