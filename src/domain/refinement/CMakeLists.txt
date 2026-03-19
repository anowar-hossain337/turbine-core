
add_library( domain-refinement )
add_library( walberla::domain-refinement ALIAS domain-refinement )
target_link_libraries( domain-refinement
        PUBLIC
        walberla::core
        walberla::blockforest
        walberla::domain_decomposition
        walberla::field
        walberla::lbm
        domain
        wind_turbine_core
)

target_sources( domain-refinement
        PRIVATE
        level_determination/CurlBasedLevelDetermination.h
        level_determination/QCriterionBasedLevelDetermination.h
        level_determination/VectorGradientBasedLevelDetermination.h
        level_determination/VorticityBasedLevelDetermination.h

        LevelDeterminationFunctors.h
        RefinementType.h
        StaticRefinementHandler.h
        TurbineRefinementHandler.h
        TurbineSelectorResetter.h
        TurbineWeightAssignmentFunctor.h
)
