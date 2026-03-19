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
//======================================================================================================================

#include "WindField.h"

#include "blockforest/StructuredBlockForest.h"
#include "core/logging/Logging.h"
#include "field/AddToStorage.h"

namespace walberla {
namespace wind {

//----------------------------------------------------------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------------------------------------------------------

WindField::WindField(const shared_ptr< StructuredBlockStorage >& blocks,
                     const WindParameters& params)
   : blocks_(blocks), params_(params)
{
   // Register the velocity field on all blocks.
   const uint_t ghostLayers = uint_t(1);
   fieldID_ = field::addToStorage< VelocityField_T >(blocks_, "windVelocity",
                                                     Vector3< real_t >(real_t(0)),
                                                     field::fzyx, ghostLayers);
}

//----------------------------------------------------------------------------------------------------------------------
// initialize
//----------------------------------------------------------------------------------------------------------------------

void WindField::initialize()
{
   WALBERLA_LOG_INFO_ON_ROOT("Initializing WindField with inflow speed "
                             << params_.inflowSpeed
                             << " in direction "
                             << params_.inflowDirection);

   for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt)
   {
      auto* field = blockIt->getData< VelocityField_T >(fieldID_);
      WALBERLA_ASSERT_NOT_NULLPTR(field);

      const auto& aabb = blockIt->getAABB();
      const real_t dz  = blocks_->dz();

      WALBERLA_FOR_ALL_CELLS_XYZ(field,
         // Compute the cell centre height in global coordinates.
         const real_t globalZ = aabb.zMin() + (real_c(z) + real_t(0.5)) * dz;

         real_t speed;
         if (globalZ > params_.roughnessLength)
            speed = params_.logLawSpeed(globalZ);
         else
            speed = real_t(0);

         field->get(x, y, z) = params_.inflowDirection * speed;
      )
   }
}

//----------------------------------------------------------------------------------------------------------------------
// update
//----------------------------------------------------------------------------------------------------------------------

void WindField::update(real_t /*dt*/)
{
   // A production implementation would apply turbulence models or LBM streaming
   // here. For the initial release we keep the field constant (steady-state
   // assumption) and rely on the LBM module for time evolution.
}

//----------------------------------------------------------------------------------------------------------------------
// getVelocity / setVelocity
//----------------------------------------------------------------------------------------------------------------------

Vector3< real_t > WindField::getVelocity(IBlock& block,
                                          cell_idx_t x,
                                          cell_idx_t y,
                                          cell_idx_t z) const
{
   const auto* field = block.getData< VelocityField_T >(fieldID_);
   WALBERLA_ASSERT_NOT_NULLPTR(field);
   return field->get(x, y, z);
}

void WindField::setVelocity(IBlock& block,
                             cell_idx_t x,
                             cell_idx_t y,
                             cell_idx_t z,
                             const Vector3< real_t >& vel)
{
   auto* field = block.getData< VelocityField_T >(fieldID_);
   WALBERLA_ASSERT_NOT_NULLPTR(field);
   field->get(x, y, z) = vel;
}

} // namespace wind
} // namespace walberla
