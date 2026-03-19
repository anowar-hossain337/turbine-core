
if( WALBERLA_BUILD_WITH_GPU_SUPPORT )

    target_link_libraries( turbine_topology
            PUBLIC
            walberla::gpu force_model force_model-gpu control_model control_model-gpu
    )

    target_sources( turbine_topology
            PRIVATE
            GPUTurbineTopology.h
            GPUTurbineTopology_GlobalFunctions.h
            Tree.h
            )

endif()
