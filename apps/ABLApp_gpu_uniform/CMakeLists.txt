waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.wd )

waLBerla_generate_target_from_python(
        NAME ABLApp_gpu_uniform_Generation
        FILE LBMGeneration.py
        OUT_FILES
        waLBerlaABL_KernelInfo.h
        waLBerlaABL_SweepCollection.cu waLBerlaABL_SweepCollection.h
        waLBerlaABL_StorageSpecification.cu waLBerlaABL_StorageSpecification.h
        waLBerlaABL_BoundaryCollection.h
        waLBerlaABL_NoSlip.cu waLBerlaABL_NoSlip.h
        waLBerlaABL_FreeSlip.cu waLBerlaABL_FreeSlip.h
        waLBerlaABL_WFB.cu waLBerlaABL_WFB.h
        waLBerlaABL_UniformUBB.cu waLBerlaABL_UniformUBB.h
        waLBerlaABL_LogLawUBB.cu waLBerlaABL_LogLawUBB.h
        waLBerlaABL_Outflow.cu waLBerlaABL_Outflow.h
        waLBerlaABL_TopOutflow.cu waLBerlaABL_TopOutflow.h
        waLBerlaABL_WelfordWFB.cu waLBerlaABL_WelfordWFB.h
        waLBerlaABL_WelfordOutput.cu waLBerlaABL_WelfordOutput.h
        waLBerlaABL_SoSResetter.cu waLBerlaABL_SoSResetter.h
        waLBerlaABL_WelfordEddyViscosity.cu waLBerlaABL_WelfordEddyViscosity.h
        waLBerlaABL_PdfSetter.cpp waLBerlaABL_PdfSetter.h
        waLBerlaABL_MacroGetter.cpp waLBerlaABL_MacroGetter.h
        waLBerlaABL_StrainRateWriter.cu waLBerlaABL_StrainRateWriter.h
        waLBerlaABL_WelfordStrainRate.cu waLBerlaABL_WelfordStrainRate.h
        FlowDriverCollection.h
        Drivers.h Drivers.cu
)

# Reason for edit start: extend the local PSM code-generation target so Phase 2 can build the real stationary-sphere coupling kernels while keeping the default ABL path available.
waLBerla_generate_target_from_python(
        NAME ABLApp_gpu_uniform_PSM_Generation
        FILE ABLPSMGeneration.py
        OUT_FILES
        waLBerlaABLPSM_Sweep.cu waLBerlaABLPSM_Sweep.h
        waLBerlaABLPSM_InitializeDomainForPSM.cu waLBerlaABLPSM_InitializeDomainForPSM.h
        waLBerlaABLPSM_KernelInfo.h
)
# Reason for edit end: extend the local PSM code-generation target so Phase 2 can build the real stationary-sphere coupling kernels while keeping the default ABL path available.

waLBerla_add_executable(
        NAME ABLApp_gpu_uniform
        FILES main.cu
        # Reason for edit start: link the shared output helper library and Phase 2 MesaPD dependencies so the app can build both generated paths while keeping runtime PSM disabled by default.
        DEPENDS walberla::core walberla::blockforest walberla::domain_decomposition walberla::field walberla::lbm
                walberla::timeloop walberla::gpu walberla::output walberla::mesa_pd walberla::lbm_mesapd_coupling
                ABLApp_gpu_uniform_Generation
                ABLApp_gpu_uniform_PSM_Generation
                domain component
                wind_turbine_core
        # Reason for edit end: link the shared output helper library and Phase 2 MesaPD dependencies so the app can build both generated paths while keeping runtime PSM disabled by default.
)
