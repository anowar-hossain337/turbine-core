
if(WALBERLA_BUILD_WITH_CUDA)

add_library( force_model-gpu )
add_library( walberla::force_model-gpu ALIAS force_model-gpu )
target_link_libraries( force_model-gpu
        PUBLIC
        force_model
        point3
        wind_turbine_core
)

target_sources( force_model-gpu
        PRIVATE
        CopyToGPU.h
)

endif()