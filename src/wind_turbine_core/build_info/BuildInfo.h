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
//! \file BuildInfo.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_BUILDINFO_H
#define TURBINECORE_BUILDINFO_H

#pragma once

namespace turbine_core {

    namespace build_info {

        const char * getGitSHA1();
        const char * getSourceDir();
        const char * getBinaryDir();

        const char * doubleAccuracy();
        const char * buildWithMPI();
        const char * buildWithOMP();
        const char * buildWithCUDA();
        const char * optimizeForLocalhost();

    } // namespace build_info


} // namespace turbine_core


#define TURBINE_CORE_GIT_SHA1   ::turbine_core::build_info::getGitSHA1()
#define TURBINE_CORE_SOURCE_DIR ::turbine_core::build_info::getSourceDir()
#define TURBINE_CORE_BUILD_DIR  ::turbine_core::build_info::getBinaryDir()


#define WALBERLA_DOUBLE_ACCURACY_STRING   ::turbine_core::build_info::doubleAccuracy()
#define WALBERLA_BUILD_WITH_MPI_STRING    ::turbine_core::build_info::buildWithMPI()
#define WALBERLA_BUILD_WITH_OPENMP_STRING ::turbine_core::build_info::buildWithOMP()
#define WALBERLA_BUILD_WITH_CUDA_STRING   ::turbine_core::build_info::buildWithCUDA()
#define WALBERLA_OPTIMIZE_FOR_LOCALHOST   ::turbine_core::build_info::optimizeForLocalhost()


#endif // TURBINECORE_BUILDINFO_H 
