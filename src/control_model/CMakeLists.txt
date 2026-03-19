
add_library( control_model )
add_library( walberla::control_model ALIAS control_model )
target_link_libraries( control_model
        PUBLIC
        walberla::core
        point3
        walberla_helper
        )

target_sources( control_model
        PRIVATE
        RAWSAnglesModel.h
        RAWSOmegaModel.h
        simpleYawController.h
        timeControlVelocities.h
        ControlModel.h
        )

add_subdirectory(gpu)
