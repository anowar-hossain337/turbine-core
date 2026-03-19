walberla_add_executable(NAME WindSimulation
                        FILES WindSimulation.cpp
                        DEPENDS blockforest core domain_decomposition field lbm
                                timeloop wind_core wind_lbm wind_boundary)
