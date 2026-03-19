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

#include "core/DataTypes.h"
#include "core/math/Vector3.h"

namespace walberla {
namespace wind {

/**
 * \brief Parameters controlling the wind simulation.
 *
 * All quantities are expressed in lattice units unless stated otherwise.
 */
struct WindParameters
{
   /// Inflow wind velocity direction (unit vector).
   Vector3< real_t > inflowDirection{ real_t(1), real_t(0), real_t(0) };

   /// Magnitude of the inflow wind speed (lattice units / time step).
   real_t inflowSpeed{ real_t(0.05) };

   /// Reference height for the wind profile (lattice units). Used in the
   /// log-law vertical wind profile.
   real_t referenceHeight{ real_t(10) };

   /// Aerodynamic roughness length z0 for the log-law wind profile.
   real_t roughnessLength{ real_t(0.1) };

   /// von Kármán constant (dimensionless, default 0.41).
   real_t vonKarmanConstant{ real_t(0.41) };

   /**
    * \brief Returns the wind speed at height \p z according to the
    *        logarithmic wind profile.
    *
    * \f$ u(z) = \frac{u_\mathrm{ref}}{\ln(z_\mathrm{ref}/z_0)}
    *            \cdot \ln\!\left(\frac{z}{z_0}\right) \f$
    *
    * \param z  Height above the ground in lattice units (must be > roughnessLength).
    * \return   Wind speed at height \p z.
    */
   real_t logLawSpeed(real_t z) const;
};

} // namespace wind
} // namespace walberla
