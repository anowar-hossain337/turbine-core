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

#include "core/DataTypes.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"

#include <memory>

namespace walberla {
namespace wind {
namespace lbm {

/**
 * \brief Couples a WindField to the waLBerla LBM solver.
 *
 * WindLBMModel provides helpers to:
 *  - initialise a PDF field from a WindField velocity profile, and
 *  - inject wind forcing into the LBM collision step via a simple
 *    body-force term.
 *
 * The template parameter \p LatticeModel_T must be a waLBerla lattice
 * model type (e.g. \c walberla::lbm::D3Q19<…>).
 */
template< typename LatticeModel_T >
class WindLBMModel
{
public:
   using PdfField_T = walberla::lbm::PdfField< LatticeModel_T >;

   /**
    * \brief Constructs the coupling between a WindField and an LBM PDF field.
    *
    * \param windField   The wind velocity field used as forcing / initialisation.
    * \param pdfFieldID  Block data ID of the LBM PDF field.
    * \param latticeModel The lattice model instance (provides relaxation rate, etc.).
    */
   WindLBMModel(const shared_ptr< WindField >& windField,
                BlockDataID pdfFieldID,
                const LatticeModel_T& latticeModel)
      : windField_(windField), pdfFieldID_(pdfFieldID), latticeModel_(latticeModel)
   {}

   /**
    * \brief Initialises the PDF field on every block using the wind velocity
    *        profile stored in the WindField.
    *
    * Populations are set to the local equilibrium \f$ f^{\rm eq} \f$ for the
    * velocity and a reference density of 1.
    *
    * \param blocks  The block storage to iterate over.
    */
   void initializePDFs(const shared_ptr< StructuredBlockStorage >& blocks);

   /**
    * \brief Applies a wind-velocity body force to the LBM PDF field.
    *
    * The forcing is implemented as an equilibrium correction: populations are
    * relaxed toward the equilibrium computed with the locally stored wind
    * velocity superimposed on the fluid velocity.
    *
    * \param block  The block on which to apply the forcing.
    */
   void applyForcing(IBlock& block);

private:
   shared_ptr< WindField > windField_;
   BlockDataID pdfFieldID_;
   LatticeModel_T latticeModel_;
};

//======================================================================================================================
// Template implementation
//======================================================================================================================

template< typename LatticeModel_T >
void WindLBMModel< LatticeModel_T >::initializePDFs(
   const shared_ptr< StructuredBlockStorage >& blocks)
{
   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      auto* pdf = blockIt->getData< PdfField_T >(pdfFieldID_);
      WALBERLA_ASSERT_NOT_NULLPTR(pdf);

      WALBERLA_FOR_ALL_CELLS_XYZ(pdf,
         const Vector3< real_t > vel = windField_->getVelocity(*blockIt, x, y, z);
         pdf->setDensityAndVelocity(x, y, z, real_t(1), vel);
      )
   }
}

template< typename LatticeModel_T >
void WindLBMModel< LatticeModel_T >::applyForcing(IBlock& block)
{
   auto* pdf = block.getData< PdfField_T >(pdfFieldID_);
   WALBERLA_ASSERT_NOT_NULLPTR(pdf);

   WALBERLA_FOR_ALL_CELLS_XYZ(pdf,
      const Vector3< real_t > windVel = windField_->getVelocity(block, x, y, z);
      Vector3< real_t > fluidVel;
      real_t rho = pdf->getDensityAndVelocity(fluidVel, x, y, z);

      // Add wind contribution scaled by the relaxation parameter.
      const real_t omega = latticeModel_.collisionModel().omega();
      const Vector3< real_t > forcedVel = fluidVel + omega * windVel;

      // Reset to equilibrium with the forced velocity (simple forcing scheme).
      pdf->setDensityAndVelocity(x, y, z, rho, forcedVel);
   )
}

} // namespace lbm
} // namespace wind
} // namespace walberla
