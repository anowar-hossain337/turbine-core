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
//! \file WindTurbineDefines.h
//! \author Ani ANCIAUX-SEDRAKIAN <ani.anciaux-sedrakian@ifpen.fr>
//! \author Frederic BLONDEL <frederic.blondel@ifpen.fr>
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#pragma once

#include <field/interpolators/all.h>
#include <field/distributors/all.h>

#include "walberla_helper/field/all.h"

#include "waLBerlaWind_KernelInfo.h"

#define WHEREAMI WALBERLA_LOG_INFO("In File " << __FILE__ << " at line " << std::to_string(__LINE__) << ".")

namespace turbine_core {

    static const uint_t fieldGhostLayers = 1;

    /// define different interpolators and distributors

    template<class Kernel_T>
    using KernelFieldInterpolator_T = projectors::KernelFieldInterpolator<field::Field<real_t>, Kernel_T>;

    template<class Kernel_T>
    using KernelDistributor_T = projectors::KernelFieldDistributor<field::Field<real_t>, Kernel_T>;

}
