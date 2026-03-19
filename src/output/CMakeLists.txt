
add_library( output )
add_library( walberla::output ALIAS output )
target_link_libraries( output
        PUBLIC
        walberla::core
        walberla::blockforest
        walberla::domain_decomposition
        walberla::field
        walberla::lbm
        walberla::timeloop
        walberla::vtk
        wind_turbine_core
)

target_sources( output
        PRIVATE
        LineOutput.h
        LineOutputLine.h
        LineOutputLineBase.h
        LineOutputMacroscopicLine.h
        Logging.h
        PlaneInclusionFilter.h
        TurbineBeforeFunctors.h
        TurbineOutput.h
)


if( WALBERLA_BUILD_WITH_CUDA )

    target_link_libraries( output
        PUBLIC
        cuda
    )

endif()