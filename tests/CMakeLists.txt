walberla_compile_test(FILES WindParametersTest.cpp DEPENDS core wind_core)
walberla_execute_test(NAME WindParametersTest)

walberla_compile_test(FILES WindFieldTest.cpp DEPENDS blockforest core domain_decomposition field wind_core)
walberla_execute_test(NAME WindFieldTest)
