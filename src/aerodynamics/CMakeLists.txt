
add_library( aerodynamics )
add_library( walberla::aerodynamics ALIAS aerodynamics )
target_link_libraries( aerodynamics
        PUBLIC
        walberla::core
        wind_turbine_core
)

target_sources( aerodynamics
        PRIVATE
        AirfoilPolar.h
)
