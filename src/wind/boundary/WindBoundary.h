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

#pragma once

#include "wind/WindField.h"

#include "boundary/Boundary.h"
#include "boundary/BoundaryHandling.h"
#include "core/DataTypes.h"

namespace walberla {
namespace wind {
namespace boundary {

/**
 * \brief Inflow boundary condition that prescribes the wind velocity profile.
 *
 * WindInflowBoundary sets populations on inflow faces to the local equilibrium
 * distribution function using the velocity obtained from the WindField.  The
 * boundary is suitable for the inlet face of a wind simulation domain (e.g.
 * the west face when wind blows in the +x direction).
 *
 * \tparam LatticeModel_T  waLBerla lattice model (e.g. D3Q19<…>).
 * \tparam FlagField_T     waLBerla flag field type used by the BoundaryHandling.
 * \tparam flag_t          Flag value type (default: uint32_t).
 */
template< typename LatticeModel_T,
          typename FlagField_T,
          typename flag_t = uint32_t >
class WindInflowBoundary
{
public:
   static const bool threadsafe = true;

   using PdfField_T = walberla::lbm::PdfField< LatticeModel_T >;

   /**
    * \param id        String identifier for this boundary condition.
    * \param windField The wind velocity field providing inflow velocities.
    * \param flagUID   Flag UID used to mark inflow boundary cells.
    */
   WindInflowBoundary(const std::string& id,
                      const shared_ptr< WindField >& windField,
                      const FlagUID& flagUID)
      : id_(id), windField_(windField), flagUID_(flagUID)
   {}

   /// Returns the string identifier.
   const std::string& id() const { return id_; }

   /// Returns the flag UID associated with this boundary.
   const FlagUID& getFlagUID() const { return flagUID_; }

   /**
    * \brief Applies the inflow boundary condition on \p block.
    *
    * For each boundary cell flagged as inflow, the populations are set to the
    * equilibrium distribution evaluated with the local wind velocity and a
    * reference density of 1.
    *
    * \param block       Block that contains the boundary cells.
    * \param pdfFieldID  Block data ID of the PDF field.
    */
   void operator()(IBlock& block, BlockDataID pdfFieldID);

private:
   std::string id_;
   shared_ptr< WindField > windField_;
   FlagUID flagUID_;
};

//======================================================================================================================
// Template implementation
//======================================================================================================================

template< typename LatticeModel_T, typename FlagField_T, typename flag_t >
void WindInflowBoundary< LatticeModel_T, FlagField_T, flag_t >::operator()(
   IBlock& block, BlockDataID pdfFieldID)
{
   auto* pdf = block.getData< PdfField_T >(pdfFieldID);
   WALBERLA_ASSERT_NOT_NULLPTR(pdf);

   WALBERLA_FOR_ALL_CELLS_XYZ(pdf,
      const Vector3< real_t > windVel = windField_->getVelocity(block, x, y, z);
      // Prescribe equilibrium populations at the inflow face.
      pdf->setDensityAndVelocity(x, y, z, real_t(1), windVel);
   )
}

} // namespace boundary
} // namespace wind
} // namespace walberla
