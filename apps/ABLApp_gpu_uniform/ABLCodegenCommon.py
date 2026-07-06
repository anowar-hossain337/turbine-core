# Reason for edit start: centralize the shared ABL code-generation settings so the plain solver and the future local PSM generator stay numerically aligned from day one.
from pathlib import Path
import re

import numpy as np
import pystencils as ps
import sympy as sp
from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil, ForceModel
from lbmpy.advanced_streaming.utility import is_inplace
from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy.turbulence_models import SubgridScaleModel
from pystencils import Target
from pystencils.fast_approximation import insert_fast_divisions, insert_fast_sqrts
from pystencils.typing import TypedSymbol

ABL_KERNEL_INFO_TEMPLATE = """

using flag_t = walberla::uint8_t;
using FlagField_T = walberla::FlagField<flag_t>;

using StorageSpecification_T = walberla::lbm::{name}_StorageSpecification;
using Stencil_T              = StorageSpecification_T::Stencil;
using CommunicationStencil_T = StorageSpecification_T::CommunicationStencil;

using PdfField_T           = walberla::lbm_generated::PdfField< StorageSpecification_T >;
using PdfGPUField_T        = walberla::lbm_generated::GPUPdfField< StorageSpecification_T >;
using BoundaryCollection_T = walberla::lbm::{name}_BoundaryCollection< FlagField_T >;

using SweepCollection_T = walberla::lbm::{name}_SweepCollection;

static const walberla::FlagUID FluidFlagUID("Fluid Flag");
static const walberla::FlagUID NoSlipFlagUID("{noslip_flag}");
static const walberla::FlagUID WFBFlagUID("{wfb_flag}");
static const walberla::FlagUID SymmetryFlagUID("{symmetry_flag}");
static const walberla::FlagUID UniformInflowFlagUID("{uniform_inflow_flag}");
static const walberla::FlagUID LogLawInflowFlagUID("{loglaw_inflow_flag}");
static const walberla::FlagUID OutflowFlagUID("{outflow_flag}");
static const walberla::FlagUID TopOutflowFlagUID("{top_outflow_flag}");

namespace codegen {{

    struct KernelInfo {{

        static constexpr char stencil[]     = "{stencil}";
        static constexpr char method[]      = "{method}";
        static constexpr char forceModel[]  = "{forceModel}";

        static constexpr walberla::uint_t q = {q};

        static constexpr walberla::field::Layout layout = walberla::field::{layout};
        static constexpr char streamingPattern[] = "{streaming_pattern}";

        static constexpr bool compressible = {compressible};
        static constexpr bool zeroCentered = {zeroCentered};
        static constexpr char subgridScaleModel[] = "{subgridScaleModel}";

        static constexpr char cpuVectoriseInfo[] = "{cpuVectoriseInfo}";
        static constexpr char lbmOptimisationDict[] = "{lbmOptimisation}";

    }};

    constexpr char KernelInfo::stencil[];
    constexpr char KernelInfo::method[];
    constexpr char KernelInfo::forceModel[];

    constexpr char KernelInfo::cpuVectoriseInfo[];
    constexpr char KernelInfo::lbmOptimisationDict[];

}} // namespace codegen

"""

COMMON_FLAG_UIDS = {
    "noslip": "NoSlip Flag",
    "wfb": "WFB Flag",
    "symmetry": "Symmetry Flag",
    "uniform_inflow": "Uniform Inflow Flag",
    "loglaw_inflow": "LogLaw Inflow Flag",
    "outflow": "Outflow Flag",
    "top_outflow": "TopOutflow Flag",
}

CPU_VECTORIZE_INFO = {"nontemporal": True}
COMPILE_TIME_BLOCK_SIZE = True
MAX_THREADS = 256
GHOST_LAYERS = 1
ABL_MAX_PARTICLES_PER_CELL = 2

if COMPILE_TIME_BLOCK_SIZE:
    SWEEP_BLOCK_SIZE = (128, 1, 1)
else:
    SWEEP_BLOCK_SIZE = (
        TypedSymbol("gpuBlockSize0", np.int32),
        TypedSymbol("gpuBlockSize1", np.int32),
        TypedSymbol("gpuBlockSize2", np.int32),
    )

GPU_INDEXING_PARAMS = {"block_size": SWEEP_BLOCK_SIZE}


def _parse_bool_token(token: str):
    value = token.strip().lower()
    if value in {"1", "true", "yes", "on"}:
        return True
    if value in {"0", "false", "no", "off"}:
        return False
    return None


def load_compressible_from_prm(default=True):
    candidates = [Path.cwd() / "input.prm", Path(__file__).with_name("input.prm")]
    for prm_path in candidates:
        if not prm_path.exists():
            continue

        content = prm_path.read_text(encoding="utf-8", errors="ignore")
        parameters_match = re.search(r"Parameters\s*\{(?P<body>.*?)\}", content, flags=re.DOTALL)
        search_space = parameters_match.group("body") if parameters_match else content
        key_match = re.search(r"\bcompressible\b\s+([^\s;]+)", search_space)
        if not key_match:
            continue

        parsed = _parse_bool_token(key_match.group(1))
        if parsed is not None:
            return parsed

    return default


def create_abl_common_setup(ctx):
    target = Target.GPU
    layout = "fzyx"
    streaming_pattern = "pull"
    omega = sp.Symbol("omega")
    compressible_flow = load_compressible_from_prm(default=True)

    data_type = "double" if ctx.double_accuracy else "float32"

    density_field = ps.fields(f"density: {data_type}[3D]", layout=layout)
    velocity_field = ps.fields(f"velocity(3): {data_type}[3D]", layout=layout)
    mean_velocity_field = ps.fields(f"mean_velocity(3): {data_type}[3D]", layout=layout)
    sum_of_squares_field = ps.fields(f"sum_of_squares(9): {data_type}[3D]", layout=layout)
    force_field = ps.fields(f"force(3): {data_type}[3D]", layout=layout)
    omega_field = ps.fields(f"omega_out: {data_type}[3D]", layout=layout)
    eddy_viscosity_field = ps.fields(f"eddy_viscosity: {data_type}[3D]", layout=layout)
    mean_eddy_viscosity_field = ps.fields(f"mean_eddy_viscosity: {data_type}[3D]", layout=layout)
    strain_rate_field = ps.fields(f"strain_rate(9): {data_type}[3D]", layout=layout)
    mean_strain_rate_field = ps.fields(f"mean_strain_rate(9): {data_type}[3D]", layout=layout)

    stencil = LBStencil(Stencil.D3Q19)
    q = stencil.Q
    use_galilean_correction = stencil.name == "D3Q27"
    fourth_order_correction = 0.1 if use_galilean_correction else 0.0

    pdfs, pdfs_tmp = ps.fields(f"pdfs({q}), pdfs_tmp({q}): {data_type}[3D]", layout=layout)
    macroscopic_fields = {"density": density_field, "velocity": velocity_field}

    lbm_config = LBMConfig(
        stencil=stencil,
        streaming_pattern=streaming_pattern,
        method=Method.CUMULANT,
        relaxation_rate=omega,
        galilean_correction=use_galilean_correction,
        fourth_order_correction=fourth_order_correction,
        force_model=ForceModel.GUO,
        force=force_field.center_vector,
        subgrid_scale_model=SubgridScaleModel.SMAGORINSKY,
        compressible=compressible_flow,
        zero_centered=compressible_flow,
        omega_output_field=omega_field,
        eddy_viscosity_field=eddy_viscosity_field,
        output=macroscopic_fields,
    )

    lbm_optimisation = LBMOptimisation(
        field_layout=layout,
        symbolic_field=pdfs,
        cse_global=True,
        cse_pdfs=False,
    )

    if not is_inplace(streaming_pattern):
        field_swaps = [(pdfs, pdfs_tmp)]
    else:
        field_swaps = []

    collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_optimisation)
    collision_rule = insert_fast_divisions(collision_rule)
    collision_rule = insert_fast_sqrts(collision_rule)

    return {
        "target": target,
        "layout": layout,
        "streaming_pattern": streaming_pattern,
        "omega": omega,
        "data_type": data_type,
        "density_field": density_field,
        "velocity_field": velocity_field,
        "mean_velocity_field": mean_velocity_field,
        "sum_of_squares_field": sum_of_squares_field,
        "force_field": force_field,
        "omega_field": omega_field,
        "eddy_viscosity_field": eddy_viscosity_field,
        "mean_eddy_viscosity_field": mean_eddy_viscosity_field,
        "strain_rate_field": strain_rate_field,
        "mean_strain_rate_field": mean_strain_rate_field,
        "stencil": stencil,
        "q": q,
        "pdfs": pdfs,
        "pdfs_tmp": pdfs_tmp,
        "macroscopic_fields": macroscopic_fields,
        "lbm_config": lbm_config,
        "lbm_optimisation": lbm_optimisation,
        "field_swaps": field_swaps,
        "collision_rule": collision_rule,
        "lb_method": collision_rule.method,
        "cpu_vectorise_info": CPU_VECTORIZE_INFO,
    }


def build_abl_info_header_params(name, common_setup):
    lbm_config = common_setup["lbm_config"]
    return {
        "name": name,
        "stencil": lbm_config.stencil.name,
        "q": common_setup["q"],
        "method": type(lbm_config.method).__name__,
        "forceModel": type(lbm_config.force_model).__name__,
        "layout": common_setup["layout"],
        "streaming_pattern": common_setup["streaming_pattern"],
        "compressible": "true" if lbm_config.compressible else "false",
        "zeroCentered": "true" if lbm_config.zero_centered else "false",
        "subgridScaleModel": lbm_config.subgrid_scale_model if lbm_config.subgrid_scale_model else "false",
        "cpuVectoriseInfo": str(common_setup["cpu_vectorise_info"]),
        "lbmOptimisation": str(vars(common_setup["lbm_optimisation"])),
        "noslip_flag": COMMON_FLAG_UIDS["noslip"],
        "wfb_flag": COMMON_FLAG_UIDS["wfb"],
        "symmetry_flag": COMMON_FLAG_UIDS["symmetry"],
        "uniform_inflow_flag": COMMON_FLAG_UIDS["uniform_inflow"],
        "loglaw_inflow_flag": COMMON_FLAG_UIDS["loglaw_inflow"],
        "outflow_flag": COMMON_FLAG_UIDS["outflow"],
        "top_outflow_flag": COMMON_FLAG_UIDS["top_outflow"],
    }


def build_abl_field_typedefs(common_setup):
    data_type = common_setup["data_type"]
    layout = common_setup["layout"]
    return {
        "ScalarField_T": ps.fields(f"dummy: {data_type}[3D]", layout=layout),
        "VectorField_T": ps.fields(f"dummy(3): {data_type}[3D]", layout=layout),
        "SecondOrderTensorField_T": ps.fields(f"dummy(9): {data_type}[3D]", layout=layout),
        "ThirdOrderTensorField_T": ps.fields(f"dummy(27): {data_type}[3D]", layout=layout),
    }


def build_abl_additional_headers():
    return {
        "field/Layout.h",
        "lbm_generated/field/PdfField.h",
        "lbm_generated/field/AddToStorage.h",
        "lbm_generated/gpu/GPUPdfField.h",
        "lbm_generated/gpu/AddToStorage.h",
        "gpu/AddGPUFieldToStorage.h",
        "gpu/HostFieldAllocator.h",
        "gpu/communication/MemcpyPackInfo.h",
        "gpu/ShiftedPeriodicity.h",
        "lbm_generated/gpu/UniformGeneratedGPUPdfPackInfo.h",
        "gpu/communication/UniformGPUScheme.h",
    }
# Reason for edit end: centralize the shared ABL code-generation settings so the plain solver and the future local PSM generator stay numerically aligned from day one.
