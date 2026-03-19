##################
### BUILD_INFO ###
##################

# Get Git hash
set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${walberla_SOURCE_DIR}/cmake )
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC TurbineCore_GIT_SHA1)

# Configure build info header
configure_file ( build_info/BuildInfo.in.cpp build_info/BuildInfo.cpp )

add_library( wind_turbine_core )
add_library( walberla::wind_turbine_core ALIAS wind_turbine_core )
target_link_libraries( wind_turbine_core
        PUBLIC
        walberla::core
)

target_sources( wind_turbine_core
        PRIVATE
        build_info/BuildInfo.cpp build_info/BuildInfo.h

        math/AABB.h
        math/AABBFwd.h
        math/GenericAABB.h
        math/Matrix2.h
        math/Matrix3.h
        math/Overloads.h
        math/Quaternion.h
        math/Round.h
        math/Vector2.h
        math/Vector3.h

        utility/CenterHelper.h
        utility/StringManipulation.h

        Enums.h
        ProjectDefines.h
        Selectors.h
        WalberlaDataTypes.h
)