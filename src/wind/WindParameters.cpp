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

#include "WindParameters.h"

#include "core/math/Math.h"

#include <cmath>

namespace walberla {
namespace wind {

real_t WindParameters::logLawSpeed(real_t z) const
{
   WALBERLA_ASSERT_GREATER(z, roughnessLength,
      "Height z must be greater than the aerodynamic roughness length z0.");

   const real_t uRef = inflowSpeed;
   const real_t logZ  = std::log(z              / roughnessLength);
   const real_t logZr = std::log(referenceHeight / roughnessLength);

   return uRef * (logZ / logZr);
}

} // namespace wind
} // namespace walberla
