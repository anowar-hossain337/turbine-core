
add_library( component )
add_library( walberla::component ALIAS component)
target_link_libraries( component
        PUBLIC
        walberla::core
        discretisation
        force_model
        mpi
        point3
        wind_turbine_core
)

target_sources( component
        PRIVATE
        Component.h
        ComponentType.h
)