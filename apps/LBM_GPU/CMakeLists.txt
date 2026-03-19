
waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( airfoil_data )
waLBerla_link_files_to_builddir( *.blade )
waLBerla_link_files_to_builddir( *.tower )
waLBerla_link_files_to_builddir( *.ctrl )
waLBerla_link_files_to_builddir( *.wd )

waLBerla_generate_target_from_python(
        NAME WindTurbineCustomLatticeModelGPU_generated
        FILE WindTurbineLBMGeneration.py
        OUT_FILES
        waLBerlaWind_KernelInfo.h
        waLBerlaWind_CollectionSweepCollection.cu waLBerlaWind_CollectionSweepCollection.h
        waLBerlaWind_CollectionStorageSpecification.cu waLBerlaWind_CollectionStorageSpecification.h
        waLBerlaWind_CollectionBoundaryCollection.h
        waLBerlaWind_NoSlip.cu waLBerlaWind_NoSlip.h
        waLBerlaWind_FreeSlip.cu waLBerlaWind_FreeSlip.h
        waLBerlaWind_WFB.cu waLBerlaWind_WFB.h
        waLBerlaWind_UniformUBB.cu waLBerlaWind_UniformUBB.h
        waLBerlaWind_LogLawUBB.cu waLBerlaWind_LogLawUBB.h
        waLBerlaWind_Outflow.cu waLBerlaWind_Outflow.h
        waLBerlaWind_Welford.cu waLBerlaWind_Welford.h
        waLBerlaWind_WelfordSOP.cu waLBerlaWind_WelfordSOP.h
        waLBerlaWind_SoSResetter.cu waLBerlaWind_SoSResetter.h
        waLBerlaWind_MacroSetter.cpp waLBerlaWind_MacroSetter.h
        waLBerlaWind_MacroGetter.cpp waLBerlaWind_MacroGetter.h
        waLBerlaWind_PdfPackInfo.cu waLBerlaWind_PdfPackInfo.h
        waLBerlaWind_EddyViscosityWriter.cu waLBerlaWind_EddyViscosityWriter.h
        waLBerlaWind_StrainRateWriter.cu waLBerlaWind_StrainRateWriter.h
        FlowDriverCollection.h
        Drivers.h Drivers.cu
)

waLBerla_add_executable(

        NAME TurbineLBM_GPU
        FILES main.cu
        DEPENDS walberla::core walberla::blockforest walberla::domain_decomposition walberla::field walberla::gpu
                walberla::lbm walberla::lbm_generated walberla::timeloop
                WindTurbineCustomLatticeModelGPU_generated
                domain component external_forces turbine_topology force_model
                wind_turbine_core walberla_helper
)
