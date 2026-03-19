
add_library( mpi )
add_library( walberla::mpi ALIAS mpi )
target_link_libraries( mpi
        PUBLIC
        walberla::blockforest
        walberla::core
        wind_turbine_core
)

target_sources( mpi
        PRIVATE
        CPUMemoryAllocator.h
        CustomMemoryBuffer.h
        SynchronisationMode.h
)

if(WALBERLA_BUILD_WITH_CUDA)

    target_link_libraries( mpi
        PUBLIC cuda
    )

    target_sources( mpi
        PRIVATE GPUMemoryAllocator.h
    )

endif()
