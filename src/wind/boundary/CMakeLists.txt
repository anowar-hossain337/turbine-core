set(WIND_BOUNDARY_HEADERS
    WindBoundary.h
)

walberla_add_module(DEPENDS blockforest core domain_decomposition field boundary wind_core)
