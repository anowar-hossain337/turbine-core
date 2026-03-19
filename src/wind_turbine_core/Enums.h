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
//! \file FieldEnum.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_ENUMS_H
#define TURBINECORE_ENUMS_H

#pragma once

#include <string>
#include <algorithm>

#include <core/Abort.h>

#include "utility/StringManipulation.h"

namespace turbine_core {

    namespace creator {

        class AeroModel {
        public:
            enum Types {
                ROTOR_DISK,
                BLADES
            };

            static Types toType(const std::string & identifier) {

                auto id = string::toLowercase(identifier);
                id = string::removeSpaces(id);

                if (id == "rotordisk" || id == "actuatordisk") {
                    return Types::ROTOR_DISK;
                } else if (id == "blades" || id == "actuatorline") {
                    return Types::BLADES;
                } else {
                    WALBERLA_ABORT("Invalid aero type (" << identifier << ").")
                }

            }

            static std::string toString(const AeroModel::Types & type) {

                switch (type) {
                    case Types::ROTOR_DISK :
                        return "rotordisk";
                    case Types::BLADES :
                        return "blades";
                }

            }

        };

    }

    namespace output {

        class Fields {

        public:

            enum Types {
                PDF,
                DENSITY,
                VELOCITY,
                MEAN_VELOCITY_WFB,
                MEAN_VELOCITY_OUTPUT,
                FORCE,
                FLAG,
                OMEGA,
                EDDY_VISCOSITY,
                MEAN_EDDY_VISCOSITY,
                STRAIN_RATE,
                MEAN_STRAIN_RATE,
                SUM_OF_SQUARES,
                SUM_OF_CUBES
            };

            static Types toType(const std::string & identifier) {

                auto id = string::toLowercase(identifier);
                id = string::removeSpaces(id);

                if (id == "pdf" || id == "pdfs" || id == "pdffield") {
                    return Types::PDF;
                } else if (id == "density") {
                    return Types::DENSITY;
                } else if (id == "velocity") {
                    return Types::VELOCITY;
                } else if (id == "meanvelocitywfb") {
                    return Types::MEAN_VELOCITY_WFB;
                } else if (id == "meanvelocity" || id == "meanvelocityoutput") {
                    return Types::MEAN_VELOCITY_OUTPUT;
                } else if (id == "force") {
                    return Types::FORCE;
                } else if (id == "flag" || id == "flags" || id == "flagfield") {
                    return Types::FLAG;
                } else if (id == "omega" || id == "relaxationrate") {
                    return Types::OMEGA;
                } else if (id == "eddyviscosity" || id == "nu_t") {
                    return Types::EDDY_VISCOSITY;
                } else if (id == "meaneddyviscosity" || id == "meannu_t") {
                    return Types::MEAN_EDDY_VISCOSITY;
                } else if (id == "strain" || id == "strainrate") {
                    return Types::STRAIN_RATE;
                } else if (id == "meanstrain" || id == "meanstrainrate") {
                    return Types::MEAN_STRAIN_RATE;
                } else if (id == "sumofsquares" || id == "m2"){
                    return Types::SUM_OF_SQUARES;
                } else if (id == "sumofcubes" || id == "m3"){
                    return Types::SUM_OF_CUBES;
                } else {
                    WALBERLA_ABORT("Invalid field name (" << identifier << ").")
                }

            }

            static std::string toString(const Fields::Types & type) {

                switch (type) {
                    case Types::PDF :
                        return "pdfField";
                    case Types::DENSITY :
                        return "density";
                    case Types::VELOCITY :
                        return "velocity";
                    case Types::MEAN_VELOCITY_WFB :
                        return "mean velocity WFB";
                    case Types::MEAN_VELOCITY_OUTPUT :
                        return "mean velocity";
                    case Types::FORCE :
                        return "force";
                    case Types::FLAG :
                        return "flag field";
                    case Types::OMEGA:
                        return "omega";
                    case Types::EDDY_VISCOSITY:
                        return "eddy viscosity";
                    case Types::MEAN_EDDY_VISCOSITY:
                        return "mean eddy viscosity";
                    case Types::STRAIN_RATE:
                        return "strain rate";
                    case Types::MEAN_STRAIN_RATE:
                        return "mean strain rate";
                    case Types::SUM_OF_SQUARES:
                        return "sum of squares";
                    case Types::SUM_OF_CUBES:
                        return "sum of cubes";
                }

            }

        };

    } // namespace output

} // namespace turbine_core

#endif // TURBINECORE_ENUMS_H
