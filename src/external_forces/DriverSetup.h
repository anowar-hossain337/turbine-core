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
//! \file DriverSetup.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef TURBINECORE_DRIVERSETUP_H
#define TURBINECORE_DRIVERSETUP_H

#include <functional>

#include <core/Set.h>
#include <core/selectable/IsSetSelected.h>

#include <core/uid/SUID.h>

#include <domain_decomposition/IBlock.h>
#include <blockforest/BlockForest.h>
#include <blockforest/StructuredBlockForest.h>

#include "wind_turbine_core/utility/StringManipulation.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

namespace walberla::domain_decomposition {
    class IBlock;
}

namespace turbine_core {

    namespace external_forces {

        class DriverSetup {

        public:

            enum Type {
                None,
                PressureGradient,
                DynamicPressureGradient,
                CoriolisForce
            };

            static Type toType(const std::string &identifier) {

                auto id = string::toLowercase(identifier);
                id = string::removeSpaces(id);

                if (id == "none") {
                    return Type::None;
                } else if (id == "pressure" || id == "pressuredriven" || id == "pressuregradient") {
                    return Type::PressureGradient;
                } else if (id == "dynamicpressure" || id == "dynamicpressuredriven" || id == "dynamicpressuregradient") {
                    return Type::DynamicPressureGradient;
                } else if (id == "coriolis" || id == "coriolisdriven" || id == "coriolisforce") {
                    return Type::CoriolisForce;
                } else {
                    WALBERLA_ABORT("Invalid driver setup name (" << identifier << ").")
                }

            }

            static std::string toString(const DriverSetup::Type &type) {

                switch (type) {
                    case Type::None :
                        return "None";
                    case Type::PressureGradient :
                        return "PressureGradient";
                    case Type::DynamicPressureGradient :
                        return "DynamicPressureGradient";
                    case Type::CoriolisForce :
                        return "CoriolisForce";
                }

            }

        };

        struct Driver {
            virtual void driveFlow(walberla::IBlock * block, const uint_t, const uint_t) = 0;
        };

    } // namespace external_forces

} // namespace turbine_core

#endif // TURBINECORE_DRIVERSETUP_H