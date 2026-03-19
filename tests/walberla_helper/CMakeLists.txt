#################
### BlockInfo ###
#################

waLBerla_compile_test( FILES blockforest/BlockInfoTest.cpp )
waLBerla_execute_test( NAME BlockInfoTest )


#############
### Field ###
#############

waLBerla_compile_test( FILES field/TurbineCoreFieldTest.cpp )
waLBerla_execute_test( NAME TurbineCoreFieldTest )

########################################
### Macroscopic Variable Calculation ###
########################################

waLBerla_compile_test( FILES field/TurbineCoreMacroscopicVariableCalculatorTest.cpp )
waLBerla_execute_test( NAME TurbineCoreMacroscopicVariableCalculatorTest )

#####################
### Interpolators ###
#####################

waLBerla_compile_test( FILES field/interpolators/TurbineCoreFieldInterpolationTest.cpp )
waLBerla_execute_test( NAME TurbineCoreFieldInterpolationTest )

waLBerla_compile_test( FILES field/interpolators/TurbineCoreMacroscopicFieldInterpolationTest.cpp )
waLBerla_execute_test( NAME TurbineCoreMacroscopicFieldInterpolationTest )

####################
### Distributors ###
####################

waLBerla_compile_test( FILES field/distributors/TurbineCoreFieldDistributorTest.cpp )
waLBerla_execute_test( NAME TurbineCoreFieldDistributorTest )