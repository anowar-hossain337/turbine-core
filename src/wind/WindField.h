//======================================================================================================================
//
//  This file is part of waLBerla-wind, a waLBerla extension for wind flow simulation.
//
//  Copyright (C) 2026 waLBerla-wind contributors
//
//  waLBerla-wind is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  waLBerla-wind is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with waLBerla-wind. If not, see <http://www.gnu.org/licenses/>.
//
//======================================================================================================================

#pragma once

#include "WindParameters.h"

#include "core/DataTypes.h"
#include "core/math/Vector3.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/GhostLayerField.h"

#include <memory>

namespace walberla {
namespace wind {

/**
 * \brief Represents a 3D wind velocity field on a structured block lattice.
 *
 * WindField stores a vector velocity field for each cell in the domain and
 * provides methods to initialize, update, and query wind velocities. It is
 * designed to be used within the waLBerla block-structured framework and can
 * be coupled with LBM solvers.
 */
class WindField
{
public:
   /// Type of the underlying ghost-layer velocity field (3 components per cell).
   using VelocityField_T = GhostLayerField< Vector3< real_t >, 1 >;

   /**
    * \brief Constructs a WindField attached to the given block storage.
    *
    * \param blocks    Shared pointer to the structured block storage.
    * \param params    Wind simulation parameters.
    */
   WindField(const shared_ptr< StructuredBlockStorage >& blocks, const WindParameters& params);

   ~WindField() = default;

   /// Initializes the wind field with a uniform inflow velocity from \p params.
   void initialize();

   /**
    * \brief Updates the wind field by one time step of length \p dt.
    *
    * This method applies a simple advection model: the field is shifted in the
    * inflow direction by one cell per \p dt * inflowSpeed units.
    *
    * \param dt  Time-step size (in simulation units).
    */
   void update(real_t dt);

   /**
    * \brief Returns the wind velocity at the block-local cell \p (x, y, z).
    *
    * \param block  The block that owns the cell.
    * \param x      Local x-index.
    * \param y      Local y-index.
    * \param z      Local z-index.
    * \return Wind velocity vector.
    */
   Vector3< real_t > getVelocity(IBlock& block, cell_idx_t x, cell_idx_t y, cell_idx_t z) const;

   /**
    * \brief Sets the wind velocity at the block-local cell \p (x, y, z).
    *
    * \param block  The block that owns the cell.
    * \param x      Local x-index.
    * \param y      Local y-index.
    * \param z      Local z-index.
    * \param vel    New wind velocity vector.
    */
   void setVelocity(IBlock& block, cell_idx_t x, cell_idx_t y, cell_idx_t z,
                    const Vector3< real_t >& vel);

   /// Returns a const reference to the wind parameters used by this field.
   const WindParameters& getParameters() const { return params_; }

   /// Returns the block data ID used to access the velocity field on each block.
   BlockDataID getFieldID() const { return fieldID_; }

private:
   shared_ptr< StructuredBlockStorage > blocks_;
   WindParameters params_;
   BlockDataID fieldID_;
};

} // namespace wind
} // namespace walberla
