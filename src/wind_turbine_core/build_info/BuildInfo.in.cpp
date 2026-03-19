//======================================================================================================================
/*!
 *  \file   BuildInfo.cpp
 *  \brief  Configured by CMake, contains information about the current build
 */
//======================================================================================================================

namespace turbine_core {

    namespace build_info {

        const char * getGitSHA1() { return "@TurbineCore_GIT_SHA1@"; }
        const char * getSourceDir() { return "@TurbineCore_SOURCE_DIR@"; }
        const char * getBinaryDir() { return "@TurbineCore_BINARY_DIR@"; }

        const char * doubleAccuracy() { return "@WALBERLA_DOUBLE_ACCURACY@"; }
        const char * buildWithMPI() { return "@WALBERLA_BUILD_WITH_MPI@"; }
        const char * buildWithOMP() { return "@WALBERLA_BUILD_WITH_OPENMP@"; }
        const char * buildWithCUDA() { return "@WALBERLA_BUILD_WITH_CUDA@"; }
        const char * optimizeForLocalhost() { return "@WALBERLA_OPTIMIZE_FOR_LOCALHOST@"; }

    } // namespace build_info

} // namespace turbine_core
