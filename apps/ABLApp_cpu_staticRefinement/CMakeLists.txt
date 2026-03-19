
waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.wd )

waLBerla_generate_target_from_python(
        NAME ABLApp_cpu_static_refinement_Generation
        FILE LBMGeneration.py
        OUT_FILES
        waLBerlaABL_KernelInfo.h
        waLBerlaABL_SweepCollection.cpp waLBerlaABL_SweepCollection.h
        waLBerlaABL_StorageSpecification.cpp waLBerlaABL_StorageSpecification.h
        waLBerlaABL_BoundaryCollection.h
        waLBerlaABL_NoSlip.cpp waLBerlaABL_NoSlip.h
        waLBerlaABL_FreeSlip.cpp waLBerlaABL_FreeSlip.h
        waLBerlaABL_WFB.cpp waLBerlaABL_WFB.h
        waLBerlaABL_UniformUBB.cpp waLBerlaABL_UniformUBB.h
        waLBerlaABL_LogLawUBB.cpp waLBerlaABL_LogLawUBB.h
        waLBerlaABL_Outflow.cpp waLBerlaABL_Outflow.h
        waLBerlaABL_WelfordWFB.cpp waLBerlaABL_WelfordWFB.h
        waLBerlaABL_WelfordOutput.cpp waLBerlaABL_WelfordOutput.h
        waLBerlaABL_SoSResetter.cpp waLBerlaABL_SoSResetter.h
        waLBerlaABL_SoCResetter.cpp waLBerlaABL_SoCResetter.h
        waLBerlaABL_WelfordEddyViscosity.cpp waLBerlaABL_WelfordEddyViscosity.h
        waLBerlaABL_PdfSetter.cpp waLBerlaABL_PdfSetter.h
        waLBerlaABL_MacroGetter.cpp waLBerlaABL_MacroGetter.h
        waLBerlaABL_StrainRateWriter.cpp waLBerlaABL_StrainRateWriter.h
        waLBerlaABL_WelfordStrainRate.cpp waLBerlaABL_WelfordStrainRate.h
        FlowDriverCollection.h
        Drivers.h Drivers.cpp
)

waLBerla_add_executable(
        NAME ABLApp_cpu_static_refinement
        FILES main.cpp
        DEPENDS walberla::core walberla::blockforest walberla::domain_decomposition walberla::field walberla::lbm
                walberla::timeloop
                ABLApp_cpu_static_refinement_Generation
                domain component
                wind_turbine_core
)
