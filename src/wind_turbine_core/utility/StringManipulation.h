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
//! \file StringManipulation.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_STRINGMANIPULATION_H
#define TURBINECORE_STRINGMANIPULATION_H

#pragma once

namespace turbine_core {

    namespace string {

        std::string toLowercase(const std::string & str) {
            auto ret = str;
            std::transform(str.begin(), str.end(), ret.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            return ret;
        }

        std::string removeSpaces(const std::string & str) {
            auto ret = str;
            ret.erase(std::remove_if(ret.begin(), ret.end(), [](unsigned char c) { return std::isspace(c); }),
                      ret.end());
            return ret;
        }

    }

} // namespace turbine_core

#endif // TURBINECORE_STRINGMANIPULATION_H 
