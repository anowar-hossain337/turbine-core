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
//! \file Selectors.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_SELECTORS_H
#define TURBINECORE_SELECTORS_H

#pragma once

#include <core/Set.h>
#include <core/uid/SUID.h>

namespace turbine_core {

    namespace sets {

        static walberla::Set<walberla::SUID> refinementSelector()                    { return walberla::Set<walberla::SUID>("refined"); }
        static walberla::Set<walberla::SUID> turbineSelector()                       { return walberla::Set<walberla::SUID>("turbine"); }
        static walberla::Set<walberla::SUID> turbineSelector(const uint_t number)    { return walberla::Set<walberla::SUID>("turbine" + std::to_string(number)); }
        static walberla::Set<walberla::SUID> empty()                                 { return walberla::Set<walberla::SUID>("empty"); }
        static walberla::Set<walberla::SUID> none()                                  { return walberla::Set<walberla::SUID>::emptySet(); }

    } //sets

} // turbine_core 

#endif //TURBINECORE_SELECTORS_H