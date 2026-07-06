from dataclasses import replace

# This is added to allow reading parameters from input.prm file to set compressibility of the flow. The function looks for a parameter named "compressible" in the input.prm file and parses its value as a boolean. If the parameter is not found or cannot be parsed, it defaults to True.
from pathlib import Path
import re
# The above imports are added to allow reading parameters from input.prm file to set compressibility of the flow. The function looks for a parameter named "compressible" in the input.prm file and parses its value as a boolean. If the parameter is not found or cannot be parsed, it defaults to True.

import pystencils as ps
import numpy as np
import sympy as sp
from lbmpy.macroscopic_value_kernels import pdf_initialization_assignments, macroscopic_values_getter
from lbmpy.flow_statistics import welford_assignments
from lbmpy.utils import second_order_moment_tensor
from pystencils import Target
from pystencils.fast_approximation import insert_fast_sqrts, insert_fast_divisions
from pystencils.typing import TypedSymbol
from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy.advanced_streaming import Timestep, get_timesteps
from lbmpy.advanced_streaming.utility import is_inplace
from lbmpy.boundaries import NoSlip, FreeSlip, UBB, ExtrapolationOutflow, WallFunctionBounce, MoninObukhovSimilarityTheory
from lbmpy_walberla.additional_data_handler import UBBAdditionalDataHandler
from lbmpy_walberla import generate_alternating_lbm_boundary
from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil, ForceModel
from lbmpy.turbulence_models import SubgridScaleModel

from pystencils_walberla import CodeGeneration, generate_info_header, generate_sweep
from lbmpy_walberla import generate_lbm_package, lbm_boundary_generator
from lbmpy_walberla.additional_data_handler import OutflowAdditionalDataHandler

from lbmpy.relaxationrates import *

from codegen_walberla_wind import generate_flow_driver_collection
# Reason for edit start: import the shared ABL codegen helper so this generator and the new local PSM scaffold stay locked to the same numerical setup.
from ABLCodegenCommon import (
    ABL_KERNEL_INFO_TEMPLATE,
    COMMON_FLAG_UIDS,
    GHOST_LAYERS,
    GPU_INDEXING_PARAMS,
    MAX_THREADS,
    build_abl_additional_headers,
    build_abl_field_typedefs,
    build_abl_info_header_params,
    create_abl_common_setup,
)
# Reason for edit end: import the shared ABL codegen helper so this generator and the new local PSM scaffold stay locked to the same numerical setup.

with CodeGeneration() as ctx:
    # Reason for edit start: source the ABL stencil, forcing, SGS model, and shared field layout from one helper so this generator remains behavior-identical while becoming reusable for the future PSM path.
    common_setup = create_abl_common_setup(ctx)
    target = common_setup["target"]
    layout = common_setup["layout"]
    streaming_pattern = common_setup["streaming_pattern"]
    omega = common_setup["omega"]
    data_type = common_setup["data_type"]
    density_field = common_setup["density_field"]
    velocity_field = common_setup["velocity_field"]
    mean_velocity_field = common_setup["mean_velocity_field"]
    sum_of_squares_field = common_setup["sum_of_squares_field"]
    force_field = common_setup["force_field"]
    omega_field = common_setup["omega_field"]
    eddy_viscosity_field = common_setup["eddy_viscosity_field"]
    mean_eddy_viscosity_field = common_setup["mean_eddy_viscosity_field"]
    strain_rate_field = common_setup["strain_rate_field"]
    mean_strain_rate_field = common_setup["mean_strain_rate_field"]
    stencil = common_setup["stencil"]
    q = common_setup["q"]
    pdfs = common_setup["pdfs"]
    macroscopic_fields = common_setup["macroscopic_fields"]
    lbm_config = common_setup["lbm_config"]
    lbm_optimisation = common_setup["lbm_optimisation"]
    collision_rule = common_setup["collision_rule"]
    lb_method = common_setup["lb_method"]
    cpu_vectorise_info = common_setup["cpu_vectorise_info"]
    # Reason for edit end: source the ABL stencil, forcing, SGS model, and shared field layout from one helper so this generator remains behavior-identical while becoming reusable for the future PSM path.

    # Welford update
    welford_wfb_update = welford_assignments(field=velocity_field, mean_field=mean_velocity_field)
    generate_sweep(ctx, "waLBerlaABL_WelfordWFB", welford_wfb_update, target=target,
                   gpu_indexing_params=GPU_INDEXING_PARAMS, max_threads=MAX_THREADS)

    # Welford update for output
    welford_output_update = welford_assignments(field=velocity_field, mean_field=mean_velocity_field,
                                                sum_of_squares_field=sum_of_squares_field)
    generate_sweep(ctx, "waLBerlaABL_WelfordOutput", welford_output_update, target=target,
                   gpu_indexing_params=GPU_INDEXING_PARAMS, max_threads=MAX_THREADS)

    @ps.kernel
    def sos_resetter():
        for d in range(stencil.D**2):
            field_access = sum_of_squares_field.center.at_index(d)
            field_access @= sp.Float(0)

    generate_sweep(ctx, "waLBerlaABL_SoSResetter", ps.AssignmentCollection(sos_resetter), target=target,
                   gpu_indexing_params=gpu_indexing_params, max_threads=max_threads)

    welford_nut_update = welford_assignments(field=eddy_viscosity_field, mean_field=mean_eddy_viscosity_field)
    generate_sweep(ctx, "waLBerlaABL_WelfordEddyViscosity", welford_nut_update, target=target,
                   gpu_indexing_params=gpu_indexing_params, max_threads=max_threads)
    welford_strain_update = welford_assignments(field=strain_rate_field, mean_field=mean_strain_rate_field)
    generate_sweep(ctx, "waLBerlaABL_WelfordStrainRate", welford_strain_update, target=target,
                   gpu_indexing_params=gpu_indexing_params, max_threads=max_threads)

    # PDF Setter -> used for initialisation before 0th timestep (CPU side)
    initial_rho = sp.Symbol('rho_0')
    pdfs_setter = pdf_initialization_assignments(lb_method=lb_method,
                                                 density=density_field.center,
                                                 velocity=velocity_field.center_vector,
                                                 pdfs=pdfs,
                                                 streaming_pattern=streaming_pattern, previous_timestep=get_timesteps(streaming_pattern)[0])

    generate_sweep(ctx, "waLBerlaABL_PdfSetter", pdfs_setter)

    # Macro getter -> used when loading snapshot (CPU side)
    getter_assignments = macroscopic_values_getter(lb_method, velocity=velocity_field.center_vector,
                                                   pdfs=pdfs.center_vector, density=density_field.center)
    generate_sweep(ctx, 'waLBerlaABL_MacroGetter', getter_assignments)

    # GENERATE BOUNDARIES
    noslip_uid = COMMON_FLAG_UIDS["noslip"]
    wfb_uid = COMMON_FLAG_UIDS["wfb"]
    symmetry_uid = COMMON_FLAG_UIDS["symmetry"]
    uniform_inflow_uid = COMMON_FLAG_UIDS["uniform_inflow"]
    loglaw_inflow_uid = COMMON_FLAG_UIDS["loglaw_inflow"]
    outflow_uid = COMMON_FLAG_UIDS["outflow"]
    top_outflow_uid = COMMON_FLAG_UIDS["top_outflow"]
    fix_density_uid = 'FixDensity Flag'

    free_slip = lbm_boundary_generator(class_name='waLBerlaABL_FreeSlip', flag_uid=symmetry_uid,
                                       boundary_object=FreeSlip(lb_method.stencil, normal_direction=(0, 0, -1)))

    no_slip = lbm_boundary_generator(class_name='waLBerlaABL_NoSlip', flag_uid=noslip_uid, boundary_object=NoSlip())
    wfb = lbm_boundary_generator(class_name='waLBerlaABL_WFB', flag_uid=wfb_uid,
                                 boundary_object=WallFunctionBounce(lb_method=lb_method, pdfs=pdfs, normal_direction=(0, 0, 1),
                                                                    wall_function_model=MoninObukhovSimilarityTheory(sp.Symbol("z0")),
                                                                    mean_velocity=mean_velocity_field,
                                                                    maronga_sampling_shift=sp.Symbol("sampling_shift"),
                                                                    data_type=data_type))

    uniform_ubb = lbm_boundary_generator(class_name='waLBerlaABL_UniformUBB', flag_uid=uniform_inflow_uid,
                                        boundary_object=UBB(velocity=[sp.Symbol("u_x"), sp.Symbol("u_y"), sp.Symbol("u_z")],
                                                            data_type=data_type))

    loglaw_ubb = lbm_boundary_generator(class_name='waLBerlaABL_LogLawUBB', flag_uid=loglaw_inflow_uid,
                                        boundary_object=UBB(lambda *args: None, dim=stencil.D, data_type=data_type))

    outflow_bc = ExtrapolationOutflow(normal_direction=(1, 0, 0),
                                      lb_method=lb_method, data_type=data_type,
                                      streaming_pattern=streaming_pattern, zeroth_timestep=Timestep.EVEN)
    outflow = lbm_boundary_generator(class_name='waLBerlaABL_Outflow', flag_uid=outflow_uid,
                                     boundary_object=outflow_bc,
                                     additional_data_handler=OutflowAdditionalDataHandler(lb_method.stencil, outflow_bc, target=target))

    top_outflow_bc = ExtrapolationOutflow(normal_direction=(0, 0, 1),
                                          lb_method=lb_method, data_type=data_type,
                                          streaming_pattern=streaming_pattern, zeroth_timestep=Timestep.EVEN)
    top_outflow = lbm_boundary_generator(class_name='waLBerlaABL_TopOutflow', flag_uid=top_outflow_uid,
                                         boundary_object=top_outflow_bc,
                                         additional_data_handler=OutflowAdditionalDataHandler(lb_method.stencil, top_outflow_bc, target=target))

    # strain rate writer
    @ps.kernel
    def strain_rate_writer():
        f_neq = sp.Matrix(pdfs.center_vector) - lb_method.get_equilibrium_terms()
        rho = lb_method.conserved_quantity_computation.density_symbol
        strain_rate_field.center_vector @= - 3 * omega_field.center / (2 * rho) * second_order_moment_tensor(f_neq, lb_method.stencil)

    strain_rate_ac = ps.AssignmentCollection(
        [lb_method.conserved_quantity_computation.equilibrium_input_equations_from_pdfs(pdfs.center_vector),
         *strain_rate_writer]
    )

    generate_sweep(ctx, "waLBerlaABL_StrainRateWriter", strain_rate_ac, target=target,
                   gpu_indexing_params=GPU_INDEXING_PARAMS, max_threads=MAX_THREADS)

    # FLOW DRIVERS
    generate_flow_driver_collection(ctx, "FlowDriverCollection", force_field=force_field, velocity_field=velocity_field,
                                    target=target, ghost_layers_to_include=GHOST_LAYERS, namespace='wind',
                                    gpu_indexing_params=GPU_INDEXING_PARAMS, max_threads=MAX_THREADS)

    name = 'waLBerlaABL'
    generate_lbm_package(ctx, name=f"{name}_",
                         collision_rule=collision_rule,
                         lbm_config=lbm_config, lbm_optimisation=lbm_optimisation,
                         nonuniform=False, boundaries=[no_slip, free_slip, wfb, uniform_ubb, loglaw_ubb, outflow, top_outflow],
                         macroscopic_fields=macroscopic_fields,
                         target=target, gpu_indexing_params=GPU_INDEXING_PARAMS, max_threads=MAX_THREADS,
                         cpu_vectorize_info=cpu_vectorise_info)
    # Reason for edit start: reuse the shared metadata builders so the plain ABL kernel header and the new PSM scaffold report the same numerical contract.
    info_header_params = build_abl_info_header_params(name, common_setup)
    field_typedefs = build_abl_field_typedefs(common_setup)
    additional_headers = build_abl_additional_headers()

    generate_info_header(ctx, 'waLBerlaABL_KernelInfo',
                         field_typedefs=field_typedefs,
                         additional_headers=additional_headers,
                         additional_code=ABL_KERNEL_INFO_TEMPLATE.format(**info_header_params)
                         )
    # Reason for edit end: reuse the shared metadata builders so the plain ABL kernel header and the new PSM scaffold report the same numerical contract.
