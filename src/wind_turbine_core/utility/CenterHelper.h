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
//! \file CenterHelper.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_CENTERHELPER_H
#define TURBINECORE_CENTERHELPER_H

#pragma once

#include <string>

namespace turbine_core {

    namespace utility {

        template<typename charT, typename traits = std::char_traits<charT> >
        class center_helper {
            std::basic_string<charT, traits> str_;
            const int width_;
        public:
            center_helper(std::basic_string<charT, traits> str, const int width) : str_(str), width_(width) {}
            template<typename a, typename b>
            friend std::basic_ostream<a, b>& operator<<(std::basic_ostream<a, b>& s, const center_helper<a, b>& c);
        };

        template<typename charT, typename traits = std::char_traits<charT> >
        center_helper<charT, traits> centered(std::basic_string<charT, traits> str, const int width = 0 ) {
            return center_helper<charT, traits>(str, width);
        }

        // redeclare for std::string directly so we can support anything that implicitly converts to std::string
        center_helper<std::string::value_type, std::string::traits_type> centered(const std::string& str, const int width = 0 ) {
            return center_helper<std::string::value_type, std::string::traits_type>(str, width);
        }

        template<typename charT, typename traits>
        std::basic_ostream<charT, traits>& operator<<(std::basic_ostream<charT, traits>& s, const center_helper<charT, traits>& c) {
            std::streamsize w = (!c.width_) ? s.width() : c.width_;
            if (w > int(c.str_.length())) {
                std::streamsize left = (w + c.str_.length()) / 2;
                s << "/// ";
                s.width(left - 4);
                s << c.str_;
                s.width(w - left);
                s << " ///";
            } else {
                s << "/// " << c.str_ << " ///";
            }
            return s;
        }

    } // namespace utility

} // namespace turbine_core

#endif // TURBINECORE_CENTERHELPER_H 
