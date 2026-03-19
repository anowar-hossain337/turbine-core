set(WIND_LBM_HEADERS
    WindLBMModel.h
)

walberla_add_module(DEPENDS blockforest core domain_decomposition field lbm wind_core)
