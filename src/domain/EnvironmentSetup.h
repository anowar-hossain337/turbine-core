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
//! \file EnvironmentSetup.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_ENVIRONMENTSETUP_H
#define TURBINECORE_ENVIRONMENTSETUP_H

#pragma once

#include <field/FlagUID.h>

#include "wind_turbine_core/utility/StringManipulation.h"

namespace turbine_core {

    class EnvironmentSetup {

    public:

        enum Type {
            Periodic,
            Tunnel,
            Open,
            FreeSlip
        };

        static Type toType(const std::string & identifier) {

            auto id = string::toLowercase(identifier);
            id = string::removeSpaces(id);

            if (id == "periodic") {
                return Type::Periodic;
            } else if (id == "tunnel") {
                return Type::Tunnel;
            } else if (id == "open") {
                return Type::Open;
            } else if (id == "freeslip") {
                return Type::FreeSlip;
            } else {
                WALBERLA_ABORT("Invalid environment setup name (" << identifier << ").")
            }

        }

        static std::string toString(const EnvironmentSetup::Type & type) {

            switch (type) {
                case Type::Periodic :
                    return "Periodic";
                case Type::Tunnel :
                    return "Tunnel";
                case Type::Open :
                    return "Open";
                case Type::FreeSlip :
                    return "FreeSlip";
            }

        }

    };


    class InflowSetup {

    public:

        enum Type {
            Periodic,
            ShiftedPeriodic,
            InflowUniform,
            InflowLogLaw
        };

        static Type toType(const std::string & identifier) {

            auto id = string::toLowercase(identifier);
            id = string::removeSpaces(id);

            if (id == "periodic") {
                return Type::Periodic;
            } else if (id == "shiftedperiodic") {
                return Type::ShiftedPeriodic;
            } else if (id == "uniforminflow") {
                return Type::InflowUniform;
            } else if (id == "loglawinflow") {
                return Type::InflowLogLaw;
            } else {
                WALBERLA_ABORT("Invalid inflow setup name (" << identifier << ").")
            }

        }

        static std::string toString(const InflowSetup::Type & type) {

            switch (type) {
                case Type::Periodic :
                    return "Periodic";
                case Type::ShiftedPeriodic :
                    return "ShiftedPeriodic";
                case Type::InflowUniform :
                    return "InflowUniform";
                case Type::InflowLogLaw :
                    return "InflowLogLaw";
            }

        }

    };


    class OutflowSetup {

    public:

        enum Type {
            Periodic,
            ShiftedPeriodic,
            Outflow
        };

        static Type toType(const std::string & identifier) {

            auto id = string::toLowercase(identifier);
            id = string::removeSpaces(id);

            if (id == "periodic") {
                return Type::Periodic;
            } else if (id == "shiftedperiodic") {
                return Type::ShiftedPeriodic;
            } else if (id == "outflow") {
                return Type::Outflow;
            } else {
                WALBERLA_ABORT("Invalid outflow setup name (" << identifier << ").")
            }

        }

        static std::string toString(const OutflowSetup::Type & type) {

            switch (type) {
                case Type::Periodic :
                    return "Periodic";
                case Type::ShiftedPeriodic :
                    return "ShiftedPeriodic";
                case Type::Outflow :
                    return "Outflow";
            }

        }

    };


    class WallSetup {

    public:

        enum Type {
            NoSlip,
            WFB
        };

        static Type toType(const std::string & identifier) {

            auto id = string::toLowercase(identifier);
            id = string::removeSpaces(id);

            if (id == "noslip") {
                return Type::NoSlip;
            } else if (id == "wfb") {
                return Type::WFB;
            } else {
                WALBERLA_ABORT("Invalid wall setup name (" << identifier << ").")
            }

        }

        static std::string toString(const WallSetup::Type & type) {

            switch (type) {
                case Type::NoSlip :
                    return "NoSlip";
                case Type::WFB :
                    return "WFB";
            }

        }

    };

} // turbine_core

#endif
