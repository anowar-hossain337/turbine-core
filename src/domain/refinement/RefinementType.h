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
//! \file RefinementType.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_REFINEMENTTYPE_H
#define TURBINECORE_REFINEMENTTYPE_H

#pragma once

namespace turbine_core {

    namespace refinement {

        enum class RefinementType : unsigned int {
            NONE              = 0u,
            TURBINE           = 1u << 0u,
            BOUNDARY_LAYER    = 1u << 1u,
            VELOCITY_GRADIENT = 1u << 2u,
            VORTICITY         = 1u << 3u,
            CURL              = 1u << 4u,
            Q_CRITERION       = 1u << 5u,
            ANY               = ~NONE
        };

        using Base_T = std::underlying_type_t <RefinementType>;

        HOST_DEVICE_PREFIX inline constexpr RefinementType operator~ (RefinementType a) {
            return static_cast<RefinementType>(~static_cast<Base_T>(a));
        }

        HOST_DEVICE_PREFIX inline constexpr RefinementType operator| (RefinementType a, RefinementType b) {
            return static_cast<RefinementType>(static_cast<Base_T>(a) | static_cast<Base_T>(b));
        }

        HOST_DEVICE_PREFIX inline constexpr RefinementType operator& (RefinementType a, RefinementType b) {
            return static_cast<RefinementType>(static_cast<Base_T>(a) & static_cast<Base_T>(b));
        }

        HOST_DEVICE_PREFIX inline constexpr RefinementType operator^ (RefinementType a, RefinementType b) {
            return static_cast<RefinementType>(static_cast<Base_T>(a) ^ static_cast<Base_T>(b));
        }

        HOST_DEVICE_PREFIX inline constexpr bool operator!(RefinementType a) {
            return (static_cast<Base_T>(a) == static_cast<Base_T>(RefinementType::NONE));
        }

        HOST_DEVICE_PREFIX inline constexpr RefinementType& operator|=(RefinementType & a, RefinementType b) {
            return a = static_cast<RefinementType>(static_cast<Base_T>(a) | static_cast<Base_T>(b));
        }

        HOST_DEVICE_PREFIX inline constexpr RefinementType& operator&=(RefinementType & a, RefinementType b) {
            return a = static_cast<RefinementType>(static_cast<Base_T>(a) & static_cast<Base_T>(b));
        }

    } // namespace refinement

    using refinement::RefinementType;

} // namespace turbine_core

#endif // TURBINECORE_REFINEMENTTYPE_H 
