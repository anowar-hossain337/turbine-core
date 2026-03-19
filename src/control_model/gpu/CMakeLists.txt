
if(WALBERLA_BUILD_WITH_CUDA)

    add_library( control_model-gpu )
    target_link_libraries( control_model-gpu
            PUBLIC
            control_model
            point3
            wind_turbine_core
            )

    target_sources( control_model-gpu
            PRIVATE
            CopyToGPU.h
            )

endif()