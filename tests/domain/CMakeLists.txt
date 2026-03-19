############################
### Boundary setup tests ###
############################

waLBerla_compile_test( FILES BoundarySetupPeriodicityTest.cpp )
waLBerla_execute_test( NAME BoundarySetupPeriodicityTest )

waLBerla_compile_test( FILES BoundarySetupFieldFillerTest.cpp )
waLBerla_execute_test( NAME BoundarySetupFieldFillerTest )

##########################
### Domain setup tests ###
##########################

waLBerla_compile_test( FILES DomainSetupDomainDecompositionTest.cpp )
waLBerla_execute_test( NAME DomainSetupDomainDecompositionTest )
