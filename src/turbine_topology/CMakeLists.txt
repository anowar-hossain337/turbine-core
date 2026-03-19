
add_library( turbine_topology )
add_library( walberla::turbine_topology ALIAS turbine_topology )
target_link_libraries( turbine_topology
        PUBLIC
        walberla::core
        component
        discretisation
        domain
        force_model
        mpi
)

add_subdirectory(cpu)
add_subdirectory(gpu)
