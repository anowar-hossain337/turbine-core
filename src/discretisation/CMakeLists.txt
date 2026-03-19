
add_library( discretisation )
add_library( walberla::discretisation ALIAS discretisation )
target_link_libraries( discretisation
        PUBLIC
        walberla::core
        force_model
        mpi
        point3
        wind_turbine_core
)
target_sources( discretisation
        PRIVATE
        Discretisation.h
        Line.h
        Disk.h

)

add_subdirectory(gpu)