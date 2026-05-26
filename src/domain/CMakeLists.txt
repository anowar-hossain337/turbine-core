
add_library( domain )
add_library( walberla::domain ALIAS domain )
target_link_libraries( domain
        PUBLIC
        walberla::blockforest
        walberla::core
        walberla::field
        walberla::geometry
        walberla::lbm
        # BEGIN EDIT: link the waLBerla mesh module so BoundarySetup can load and voxelize STL/OBJ obstacles.
        walberla::mesh
        # END EDIT: link the waLBerla mesh module so BoundarySetup can load and voxelize STL/OBJ obstacles.
        walberla::vtk
        domain-refinement
        wind_turbine_core
)

target_sources( domain
        PRIVATE
        BoundaryHandling.h
        BoundarySetter.h
        BoundarySetup.h
        DomainInitialiser.h
        DomainSetup.h
        EnvironmentSetup.h
        TurbineDomain.h
)

add_subdirectory( refinement )
