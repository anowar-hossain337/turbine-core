# Reason for edit start: generate an ABL-compatible Phase 2 PSM sweep and one-time PDF reinitializer so the optional MesaPD path can reuse the current ABL numerics instead of switching to the stock PSM solver setup.
from dataclasses import replace

import pystencils as ps
from lbmpy import LBMOptimisation
from lbmpy.creationfunctions import create_psm_update_rule
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.partially_saturated_cells import PSMConfig
from pystencils_walberla import CodeGeneration, generate_info_header, generate_sweep

from ABLCodegenCommon import (
    ABL_MAX_PARTICLES_PER_CELL,
    GPU_INDEXING_PARAMS,
    MAX_THREADS,
    build_abl_additional_headers,
    build_abl_field_typedefs,
    build_abl_info_header_params,
    create_abl_common_setup,
)

PSM_INFO_HEADER = """
namespace codegen {{
namespace psm {{

struct KernelInfo {{
    static constexpr char stencil[] = "{stencil}";
    static constexpr char method[] = "{method}";
    static constexpr char forceModel[] = "{forceModel}";
    static constexpr walberla::uint_t q = {q};
    static constexpr char layout[] = "{layout}";
    static constexpr char streamingPattern[] = "{streaming_pattern}";
    static constexpr bool compressible = {compressible};
    static constexpr bool zeroCentered = {zeroCentered};
    static constexpr char subgridScaleModel[] = "{subgridScaleModel}";
    static constexpr char cpuVectoriseInfo[] = "{cpuVectoriseInfo}";
    static constexpr char lbmOptimisationDict[] = "{lbmOptimisation}";
    static constexpr walberla::uint_t maxParticlesPerCell = {max_particles_per_cell};
    static constexpr bool runtimeSweepAvailable = true;
}};

constexpr char KernelInfo::stencil[];
constexpr char KernelInfo::method[];
constexpr char KernelInfo::forceModel[];
constexpr char KernelInfo::layout[];
constexpr char KernelInfo::streamingPattern[];
constexpr char KernelInfo::subgridScaleModel[];
constexpr char KernelInfo::cpuVectoriseInfo[];
constexpr char KernelInfo::lbmOptimisationDict[];

}} // namespace psm
}} // namespace codegen
"""

with CodeGeneration() as ctx:
    common_setup = create_abl_common_setup(ctx)

    target = common_setup["target"]
    data_type = common_setup["data_type"]
    layout = common_setup["layout"]
    stencil = common_setup["stencil"]
    density_field = common_setup["density_field"]
    velocity_field = common_setup["velocity_field"]
    force_field = common_setup["force_field"]
    omega_field = common_setup["omega_field"]
    eddy_viscosity_field = common_setup["eddy_viscosity_field"]
    pdfs = common_setup["pdfs"]
    pdfs_tmp = common_setup["pdfs_tmp"]
    lb_method = common_setup["lb_method"]

    particle_velocities = ps.fields(
        f"particle_v({ABL_MAX_PARTICLES_PER_CELL * stencil.D}): {data_type}[3D]",
        layout=layout,
    )
    particle_forces = ps.fields(
        f"particle_f({ABL_MAX_PARTICLES_PER_CELL * stencil.D}): {data_type}[3D]",
        layout=layout,
    )
    Bs = ps.fields(f"Bs({ABL_MAX_PARTICLES_PER_CELL}): {data_type}[3D]", layout=layout)
    B = ps.fields(f"b: {data_type}[3D]", layout=layout)

    psm_config = PSMConfig(
        fraction_field=B,
        object_velocity_field=particle_velocities,
        SC=1,
        MaxParticlesPerCell=ABL_MAX_PARTICLES_PER_CELL,
        individual_fraction_field=Bs,
        particle_force_field=particle_forces,
    )

    psm_lbm_config = replace(common_setup["lbm_config"], psm_config=psm_config)
    psm_lbm_optimisation = LBMOptimisation(
        field_layout=layout,
        symbolic_field=pdfs,
        symbolic_temporary_field=pdfs_tmp,
        cse_global=common_setup["lbm_optimisation"].cse_global,
        cse_pdfs=common_setup["lbm_optimisation"].cse_pdfs,
    )

    psm_rule = create_psm_update_rule(psm_lbm_config, psm_lbm_optimisation)
    generate_sweep(
        ctx,
        "waLBerlaABLPSM_Sweep",
        psm_rule,
        field_swaps=[(pdfs, pdfs_tmp)],
        target=target,
        gpu_indexing_params=GPU_INDEXING_PARAMS,
        max_threads=MAX_THREADS,
    )

    blended_velocity = []
    for direction in range(stencil.D):
        blended_component = velocity_field.center_vector[direction] * (1.0 - B.center)
        for particle_index in range(ABL_MAX_PARTICLES_PER_CELL):
            blended_component += (
                particle_velocities.center_vector[particle_index * stencil.D + direction]
                * Bs.center_vector[particle_index]
            )
        blended_velocity.append(blended_component)

    pdf_reinitializer = macroscopic_values_setter(
        lb_method,
        density_field.center,
        blended_velocity,
        pdfs.center_vector,
    )
    generate_sweep(
        ctx,
        "waLBerlaABLPSM_InitializeDomainForPSM",
        pdf_reinitializer,
        target=target,
        gpu_indexing_params=GPU_INDEXING_PARAMS,
        max_threads=MAX_THREADS,
    )

    info_header_params = build_abl_info_header_params("waLBerlaABLPSM", common_setup)
    info_header_params["max_particles_per_cell"] = ABL_MAX_PARTICLES_PER_CELL
    generate_info_header(
        ctx,
        "waLBerlaABLPSM_KernelInfo",
        field_typedefs=build_abl_field_typedefs(common_setup),
        additional_headers=build_abl_additional_headers(),
        additional_code=PSM_INFO_HEADER.format(**info_header_params),
    )
# Reason for edit end: generate an ABL-compatible Phase 2 PSM sweep and one-time PDF reinitializer so the optional MesaPD path can reuse the current ABL numerics instead of switching to the stock PSM solver setup.
