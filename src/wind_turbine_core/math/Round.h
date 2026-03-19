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
//! \file Round.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_ROUND_H
#define TURBINECORE_ROUND_H

#pragma once

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/ProjectDefines.h"

#include <cmath>

namespace turbine_core {

    namespace math {

        HOST_DEVICE_PREFIX inline real_t round( real_t number ) {
            return ::round(number);
        }

        template< uint_t SMARTTHRESHOLD = 0 >
        HOST_DEVICE_PREFIX real_t round( const real_t number, real_t threshold ) {

            double intpart;
            auto fracpart = real_t(modf(double(number), &intpart));

            real_t sign = (number < real_t(0.0)) ? real_t(-1.0) : real_t(1.0);

            bool roundDown;
            if (SMARTTHRESHOLD == 0) {
                roundDown = ::fabs(fracpart) <= threshold;
            } else {
                constexpr real_t shift{SMARTTHRESHOLD};
                roundDown = ::fabs(fracpart) <= ( real_t(0.5) + real_t(0.5) * (number - sign * shift) / ( real_t(1.0) + ::fabs(number - sign * shift) ) );
            }

            return roundDown ? intpart : intpart + sign;

        }

        template< uint_t MULTIPLE = 2, uint_t SMARTTHRESHOLD = 0 >
        HOST_DEVICE_PREFIX inline real_t roundToMultiple( const real_t number, const real_t threshold = real_t(0.5) ) {

            constexpr auto mult = double(MULTIPLE);
            real_t sign = (number < real_t(0.0)) ? real_t(-1.0) : real_t(1.0);

            auto fracpart = real_t(fmod(double(number), mult));
            real_t intpart = number - fracpart;

            bool roundDown;
            if (SMARTTHRESHOLD == 0) {
                roundDown = ::fabs(fracpart) <= mult * threshold;
            } else {
                constexpr real_t shift{SMARTTHRESHOLD};
                const real_t posNumberShifted = sign * number - shift;
                roundDown = ::fabs(fracpart) <= mult * ( real_t(0.5) + real_t(0.5) * posNumberShifted / ( real_t(1.0) + ::fabs(posNumberShifted) ) );
            }

            return roundDown ? intpart : (intpart + sign * mult);

        }


    } // math

} // turbine_core 

#endif //TURBINECORE_ROUND_H