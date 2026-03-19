add_subdirectory(lbm)
add_subdirectory(boundary)

set(WIND_CORE_SOURCES
    WindField.cpp
    WindParameters.cpp
)

set(WIND_CORE_HEADERS
    WindField.h
    WindParameters.h
)

walberla_add_module(DEPENDS blockforest core domain_decomposition field timeloop)
