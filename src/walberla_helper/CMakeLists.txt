
add_library( walberla_helper )
add_library( walberla::walberla_helper ALIAS walberla_helper )

target_link_libraries( walberla_helper
        PUBLIC
        walberla::core
        walberla::blockforest
        walberla::domain_decomposition
        walberla::field
        wind_turbine_core
)

target_sources( walberla_helper
        PRIVATE

        blockforest/BlockInfo.h

        field/all.h

        field/distributors/KernelFieldDistributor.h
        field/distributors/NearestNeighbourFieldDistributor.h

        field/interpolators/KernelFieldInterpolator.h
        field/interpolators/NearestNeighbourFieldInterpolator.h
        field/interpolators/TrilinearFieldInterpolator.h

        field/projector_kernels/GaussianKernel.h
        field/projector_kernels/M4ProjectionKernel.h
        field/projector_kernels/SmoothedDiracDeltaKernel.h

        field/Field.h
        field/MacroscopicVariableCalculator.h
        field/Projectors.h

        refinement/CustomRecursiveTimeStep.h
)

if( WALBERLA_BUILD_WITH_CUDA OR WALBERLA_BUILD_WITH_GPU_SUPPORT )
    target_link_libraries( walberla_helper
            PUBLIC
            walberla::gpu
    )

    target_sources( walberla_helper
            PRIVATE
            field/FieldCopyGPU.h
    )
endif()