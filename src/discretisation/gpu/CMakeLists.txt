
if(WALBERLA_BUILD_WITH_CUDA)

target_sources( discretisation
        PRIVATE
        CopyToGPU.h
)

endif()