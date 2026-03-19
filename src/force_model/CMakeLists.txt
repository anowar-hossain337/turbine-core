
add_library( force_model )
add_library( walberla::force_model ALIAS force_model )
target_link_libraries( force_model
        PUBLIC
        walberla::core
        point3
        walberla_helper
)

target_sources( force_model
        PRIVATE
        ActuatorDiskModel.h
        ActuatorLineModel.h
        ForceModel.h
        TipLossModel.h
)

add_subdirectory(gpu)
