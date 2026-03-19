
add_library( external_forces )
add_library( walberla::external_forces ALIAS external_forces )
target_link_libraries( external_forces
        PUBLIC
        walberla::core
        point3
        walberla_helper
)

target_sources( external_forces
        PRIVATE
        DriverSetup.h
        FlowDriver.h
        ZeroForce.h
        PressureGradient.h
        DynamicPressureGradient.h
        CoriolisForce.h
)

if( WALBERLA_BUILD_WITH_CUDA OR WALBERLA_BUILD_WITH_GPU_SUPPORT )
    target_link_libraries( external_forces
            PUBLIC
            walberla::gpu
    )

    target_sources( external_forces
            PRIVATE
            FlowDriverGPU.h
            ZeroForceGPU.h
            PressureGradientGPU.h
            DynamicPressureGradientGPU.h
            CoriolisForceGPU.h
    )
endif()