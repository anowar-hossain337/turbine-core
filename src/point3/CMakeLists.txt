
add_library( point3 )
add_library( walberla::point3 ALIAS point3 )
target_link_libraries( point3
        PUBLIC
        aerodynamics
        walberla::core
        wind_turbine_core
)

target_sources( point3
        PRIVATE
        Point3.cpp Point3.h
)