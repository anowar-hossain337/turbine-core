from .enums import *
from pystencils import Target
from lbmpy import Stencil
from .timeloop import *
from .welford import *

field_creation_string_snapshot = """
auto & ssh = domainSetup.snapshotHandler_;

auto vectorFieldDataHandling = std::make_shared<field::FlattenedFieldHandling<{vector_field_type}>> (blocks, fieldGhostLayers, layout, allocator);
const auto forceFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->forceFile).c_str());

{force_field_variable} = blocks->loadBlockData(forceFileName, vectorFieldDataHandling, "{force_field_name}");

auto pdfFieldDataHandling = std::make_shared< walberla::lbm_generated::internal::PdfFieldHandling< StorageSpecification_T > >( blocks, storageSpec, fieldGhostLayers, layout, allocator );

const auto pdfFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->pdfFile).c_str());
{pdf_field_variable} = blocks->loadBlockData(pdfFileName, pdfFieldDataHandling, "{pdf_field_name}");
"""

field_initialiser_string_no_snapshot = """
    std::unique_ptr<domain::DomainInitialisation> initialiser{{}};
    if (initType == domain::DomainInitialisation::UNIFORM) {{
        initialiser = std::make_unique<domain::UniformInitialisation>(blocks, initialVelocity);
    }} else if (initType == domain::DomainInitialisation::ASMUTH) {{
        initialiser = std::make_unique<domain::AsmuthInitialisation>(blocks, uTau, roughnessLength, kappa, B,
                                                                     viscosity_LU, domainSetup.domainSize_, true);
    }} else if (initType == domain::DomainInitialisation::LOG_LAW) {{
        // asmuth without perturbation velocities
        initialiser = std::make_unique<domain::AsmuthInitialisation>(blocks, uTau, roughnessLength, kappa, B,
                                                                     viscosity_LU, domainSetup.domainSize_, false);
    }}
    using PdfSetter_T = walberla::pystencils::waLBerlaWind_PdfSetter;
    PdfSetter_T setter({density_field}, {force_field}, {pdf_field}, {velocity_field});
    initialiser->setViaVelocityField<{velocity_type}, PdfSetter_T>({velocity_field}, setter);
"""


def generate_field_creation(field_map: set[Fields], target: Target, snapshot: bool, generate_wfb: bool, vtk_output: bool):

    add_to_storage_string = ("{variable} = walberla::field::addToStorage<{field_type}>("
                             "blocks, \"{name}\", real_t({default_value}), layout, fieldGhostLayers, allocator);")

    result = ["std::map<output::Fields::Types, BlockDataID> fieldMap {};"]

    if generate_wfb or vtk_output:
        result.append("uint_t initialWelfordCounter = 0;\n")

    if target == Target.CPU:
        result.append("auto allocator = std::make_shared< walberla::field::AllocateAligned<real_t, 64> >();")
    else:
        result.append("auto allocator = std::make_shared< walberla::gpu::HostFieldAllocator<real_t> >(); // use pinned memory allocator for faster CPU-GPU memory transfers")

    # order set to have consistent output
    field_list = list(field_map)
    field_list = sorted(field_list, key=lambda field: field.name)

    if snapshot:
        possibly_loaded_fields = [Fields.FORCE, Fields.PDF,
                                  Fields.MEAN_VELOCITY_WFB, Fields.MEAN_VELOCITY_OUTPUT,
                                  Fields.SUM_OF_SQUARES, Fields.SUM_OF_CUBES]
        created_fields = set[Fields]()
        loaded_fields = set[Fields]()

        # BlockDataID declaration or safe for later
        for field in field_list:
            if field not in possibly_loaded_fields:
                created_fields.add(field)
            else:
                loaded_fields.add(field)
                result.append(f"BlockDataID {get_field_variable(field, Target.CPU)};")

        created_fields = list(created_fields)
        created_fields = sorted(created_fields, key=lambda field: field.name)
        loaded_fields = list(loaded_fields)
        loaded_fields = sorted(loaded_fields, key=lambda field: field.name)

        result.append("\nif(domainSetup.snapshotHandler_->loadSnapshot) {")
        result.append(field_creation_string_snapshot.format(vector_field_type=get_field_type(Fields.FORCE, None),
                                                            force_field_variable=get_field_variable(Fields.FORCE, Target.CPU),
                                                            force_field_name=get_field_name(Fields.FORCE, Target.CPU),
                                                            pdf_field_variable=get_field_variable(Fields.PDF, Target.CPU),
                                                            pdf_field_name=get_field_name(Fields.PDF, Target.CPU)).replace("\n", "\n\t"))

        load_mean_field = Fields.MEAN_VELOCITY_WFB in field_map or Fields.MEAN_VELOCITY_OUTPUT in field_map

        if load_mean_field:
            result.append("\n\tif(!ssh->meanVelFile.empty()) {")
            result.append("\t\tconst auto meanVelFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->meanVelFile).c_str());")
            result.append("\t\tinitialWelfordCounter = ssh->welfordCounter;\n")

            if generate_wfb and Fields.MEAN_VELOCITY_WFB in field_map:
                result.append(f"\t\t{get_field_variable(Fields.MEAN_VELOCITY_WFB, Target.CPU)}    = blocks->loadBlockData(meanVelFileName, vectorFieldDataHandling, "
                              f"\"{get_field_name(Fields.MEAN_VELOCITY_WFB, Target.CPU)}\");")

            if vtk_output and Fields.MEAN_VELOCITY_OUTPUT in field_map:
                result.append(f"\t\t{get_field_variable(Fields.MEAN_VELOCITY_OUTPUT, Target.CPU)} = blocks->loadBlockData(meanVelFileName, vectorFieldDataHandling, "
                              f"\"{get_field_name(Fields.MEAN_VELOCITY_OUTPUT, Target.CPU)}\");")

                if Fields.SUM_OF_SQUARES in field_map:
                    result.append(f"\n\t\tif(!ssh->sosFile.empty()) {{")
                    result.append(f"\t\t\tauto secondOrderTensorFieldDataHandling = std::make_shared<field::FlattenedFieldHandling<{get_field_type(Fields.SUM_OF_SQUARES, None)}>> (blocks, fieldGhostLayers, layout, allocator);")
                    result.append(f"\t\t\tconst auto sosFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->sosFile).c_str());")
                    result.append(f"\t\t\t{get_field_variable(Fields.SUM_OF_SQUARES, Target.CPU)} = blocks->loadBlockData(sosFileName, secondOrderTensorFieldDataHandling, \"{get_field_name(Fields.SUM_OF_SQUARES)}\");")
                    result.append(f"\t\t}}")

                if Fields.SUM_OF_CUBES in field_map:
                    result.append(f"\t\tif(!ssh->socFile.empty()) {{")
                    result.append(f"\t\t\tauto thirdOrderTensorFieldDataHandling = std::make_shared<field::FlattenedFieldHandling<{get_field_type(Fields.SUM_OF_CUBES, None)}>> (blocks, fieldGhostLayers, layout, allocator);")
                    result.append(f"\t\t\tconst auto socFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->socFile).c_str());")
                    result.append(f"\t\t\t{get_field_variable(Fields.SUM_OF_CUBES, Target.CPU)} = blocks->loadBlockData(socFileName, thirdOrderTensorFieldDataHandling, \"{get_field_name(Fields.SUM_OF_CUBES, Target.CPU)}\");")
                    result.append(f"\t\t}}")

            result.append("\t}")

        result.append("} else {")
        for field in loaded_fields:
            default_value = 1 if field == Fields.DENSITY else 0
            if field == Fields.PDF:
                result.append(f"\t{get_field_variable(field, Target.CPU)} = "
                              f"walberla::lbm_generated::addPdfFieldToStorage(blocks, \"{get_field_name(field, Target.CPU)}\", "
                              f"storageSpec, fieldGhostLayers, layout, "
                              f"walberla::Set<walberla::SUID>::emptySet(), walberla::Set<walberla::SUID>::emptySet(), allocator);")
            else:
                result.append("\t" + add_to_storage_string.format(variable=get_field_variable(field, Target.CPU),
                                                                  field_type=get_field_type(field, None),
                                                                  name=get_field_name(field, Target.CPU),
                                                                  default_value=default_value))

        result.append("}\n")

        for field in created_fields:
            default_value = 1 if field == Fields.DENSITY else 0
            if field == Fields.FLAG:
                result.append(f"BlockDataID {get_field_variable(field, None)} = "
                              f"walberla::field::addFlagFieldToStorage<{get_field_type(Fields.FLAG)}>(blocks, "
                              f"\"{get_field_name(field, None)}\", fieldGhostLayers);")

            else:
                result.append("BlockDataID " + add_to_storage_string.format(variable=get_field_variable(field, Target.CPU),
                                                                            field_type=get_field_type(field, None),
                                                                            name=get_field_name(field, Target.CPU),
                                                                            default_value=default_value))

    else:
        for field in field_list:
            default_value = 1 if field == Fields.DENSITY else 0
            if field == Fields.PDF:
                result.append(f"BlockDataID {get_field_variable(field, Target.CPU)} = "
                              f"walberla::lbm_generated::addPdfFieldToStorage(blocks, \"{get_field_name(field, Target.CPU)}\", "
                              f"storageSpec, fieldGhostLayers, layout, "
                              f"walberla::Set<walberla::SUID>::emptySet(), walberla::Set<walberla::SUID>::emptySet(), allocator);")
            elif field == Fields.FLAG:
                result.append(f"BlockDataID {get_field_variable(field)} = "
                              f"walberla::field::addFlagFieldToStorage<{get_field_type(Fields.FLAG, None)}>(blocks, "
                              f"\"{get_field_name(field)}\", fieldGhostLayers);")

            else:
                result.append("BlockDataID " + add_to_storage_string.format(variable=get_field_variable(field, Target.CPU),
                                                                            field_type=get_field_type(field, None),
                                                                            name=get_field_name(field, Target.CPU),
                                                                            default_value=default_value))

    # add to fieldMap
    result.append("\n")
    for field in field_list:
        if field == Fields.FLAG:
            result.append(f"fieldMap[output::Fields::{field.name.upper()}] = {get_field_variable(field, None)};")
        else:
            result.append(f"fieldMap[output::Fields::{field.name.upper()}] = {get_field_variable(field, Target.CPU)};")

    return "\n".join(result)


def generate_gpu_field_creation(field_map: set[Fields]):
    generation_string = ("BlockDataID {gpu_variable} = walberla::gpu::addGPUFieldToStorage<{field_type}>( blocks, "
                         "{cpu_variable}, \"{name}\", true );")

    result = ["WALBERLA_LOG_INFO_ON_ROOT(\"Add GPU fields...\")\n"]

    # order set to have consistent output
    field_list = list(field_map)
    field_list = sorted(field_list, key=lambda field: field.name)

    for field in field_list:
        name = get_field_name(field, Target.GPU)
        cpu_var = get_field_variable(field, Target.CPU)
        gpu_var = get_field_variable(field, Target.GPU)
        type_str = get_field_type(field)
        if field == Fields.PDF:
            result.append(f"BlockDataID {gpu_var} = walberla::lbm_generated::"
                       f"addGPUPdfFieldToStorage<{type_str}, StorageSpecification_T >(blocks, "
                       f"{cpu_var}, storageSpec, \"{name}\");")
        elif field == Fields.FLAG:
            pass
        else:
            result.append(generation_string.format(gpu_variable=gpu_var, cpu_variable=cpu_var,
                                                   field_type=type_str, name=name))

    return "\n".join(result)


def generate_sweep_collection(target: Target, field_map: set[Fields], runtime_gpu_indexing: bool):
    result = []

    indent = ' ' * len("SweepCollection_T sweepCollection( ")

    result.append(f"SweepCollection_T sweepCollection( blocks, {get_field_variable(Fields.DENSITY, target)}, "
                  f"{(get_field_variable(Fields.EDDY_VISCOSITY, target)+', ' if Fields.EDDY_VISCOSITY in field_map else '')}"
                  f"{get_field_variable(Fields.FORCE, target)}, "
                  f"{(get_field_variable(Fields.OMEGA, target)+', ' if Fields.OMEGA in field_map else '')}")
    result.append(f"{indent}{get_field_variable(Fields.PDF, target)}, {get_field_variable(Fields.VELOCITY, target)}, "
                  f"{'gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2], ' if runtime_gpu_indexing else ''}omega );\n")

    result.append("// for (auto& block : *blocks) {")
    result.append("//     sweepCollection.initialise(&block);")
    result.append("//     sweepCollection.initialise(&block, fieldGhostLayers);")
    result.append("// }")
    result.append("WALBERLA_MPI_BARRIER()")
    result.append("WALBERLA_LOG_INFO_ON_ROOT(\"Initialisation done\")")

    return "\n".join(result)


def generate_field_initialisation(field_map: set[Fields], with_snapshot: bool, vtk_output: bool, generate_wfb: bool,
                                  streaming_pattern: StreamingPattern):

    result = []

    result.append("\t" + field_initialiser_string_no_snapshot.format(density_field=get_field_variable(Fields.DENSITY, Target.CPU),
                                                                     force_field=get_field_variable(Fields.FORCE, Target.CPU),
                                                                     pdf_field=get_field_variable(Fields.PDF, Target.CPU),
                                                                     velocity_field=get_field_variable(Fields.VELOCITY, Target.CPU),
                                                                     velocity_type=get_field_type(Fields.VELOCITY)))

    if vtk_output and Fields.MEAN_VELOCITY_OUTPUT in field_map:
        result.append(f"\tinitialiser->setViaVelocityField<{get_field_type(Fields.MEAN_VELOCITY_OUTPUT)}, PdfSetter_T>({get_field_variable(Fields.MEAN_VELOCITY_OUTPUT, Target.CPU)}, setter);")
    if generate_wfb and Fields.MEAN_VELOCITY_WFB in field_map:
        result.append(f"\tinitialiser->setViaVelocityField<{get_field_type(Fields.MEAN_VELOCITY_WFB)}, PdfSetter_T>({get_field_variable(Fields.MEAN_VELOCITY_WFB, Target.CPU)}, setter);")

    if with_snapshot:
        no_snapshot_init = "\n".join(result)
        result = []
        result.append("if(!domainSetup.snapshotHandler_->loadSnapshot) {\n\t\t" + no_snapshot_init + "\n} else {\n")

        result.append("\tif(domainSetup.snapshotHandler_->meanVelFile.empty()) {")
        result.append("\t\tWALBERLA_ABORT(\"Mean velocity not loaded from snapshot but also not initialised.\")")
        result.append("\t}\n")
        result.append("\t// initialise velocity field from pfd field if snapshot")

        getter_sweep = """    for( auto & block : *blocks ) {{
                    {getter_name}(&block);
                }}"""

        getter_args = (f"{get_field_variable(Fields.DENSITY, Target.CPU)}, "
                       f"{get_field_variable(Fields.FORCE, Target.CPU)}, {get_field_variable(Fields.PDF, Target.CPU)}, "
                       f"{get_field_variable(Fields.VELOCITY, Target.CPU)}")

        if not streaming_pattern.is_inplace():
            result.append(f"\twalberla::pystencils::waLBerlaWind_MacroGetter getter({getter_args});")
            result.append("\t" + getter_sweep.format(getter_name="getter"))
        else:
            result.append(f"\twalberla::pystencils::waLBerlaWind_MacroGetter_EVEN getterEven({getter_args});")
            result.append(f"\twalberla::pystencils::waLBerlaWind_MacroGetter_ODD getterOdd({getter_args});")

            result.append("\tif(domainSetup.snapshotHandler_->startingTimestep % 2 == 0) {")
            result.append("\t" + getter_sweep.format(getter_name="getterOdd"))
            result.append("\t} else {")
            result.append("\t" + getter_sweep.format(getter_name="getterEven"))
            result.append("\t}")

        result.append("}")

    return "\n".join(result)


def generate_boundary_collection(target, generate_wfb):

    indent = "\t"

    result = []
    result.append(f"const auto& inflowVelocity = boundarySetup.inflowVelocity();")
    result.append(f"auto velocityInit = boundary::velocityInit(roughnessLength, kappa, uTau);")
    result.append(f"BoundaryCollection_T boundaryCollection( ")
    result.append(f"{indent}blocks, {get_field_variable(Fields.FLAG)}, {get_field_variable(Fields.PDF, target)}, "
                  f"FluidFlagUID, {(get_field_variable(Fields.FORCE, target)+', ' if generate_wfb else '')}")

    if generate_wfb:
        result.append(f"{indent}{get_field_variable(Fields.MEAN_VELOCITY_WFB, target)}, samplingHeight,")
        result.append(f"{indent}roughnessLengthRatio * referenceHeight,")

    result.append(f"{indent}inflowVelocity[0], inflowVelocity[1], inflowVelocity[2], velocityInit" + ("," if target == Target.GPU else ''))

    if target == Target.GPU:
        result.append(f"{indent}{get_field_variable(Fields.PDF, Target.CPU)}")

    result.append(");")

    return "\n".join(result)


def generate_set_field_ids(target: Target):
    return (f"farm.setFieldIDs({get_field_variable(Fields.DENSITY, target)}, "
            f"{get_field_variable(Fields.VELOCITY, target)}, {get_field_variable(Fields.FORCE, target)});")


def generate_field_resetter(target: Target, welford_output: Welford, runtime_gpu_indexing):

    result = []

    namespace = "gpu" if target == Target.GPU else "field"
    basename = welford_output.basename[0].lower() + welford_output.basename[1:]

    gpu_indexing = ', gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]' if runtime_gpu_indexing else ''

    if welford_output.active:
        if len(welford_output.welford_fields) >= 3:
            result.append(f"\nwalberla::pystencils::waLBerlaWind_SoSResetter {basename}SosResetter({get_field_variable(welford_output.sos_field, target)}{gpu_indexing});")
        if len(welford_output.welford_fields) >= 4:
            result.append(f"\nwalberla::pystencils::waLBerlaWind_SoCResetter {basename}SocResetter({get_field_variable(welford_output.soc_field, target)}{gpu_indexing});")

    return "\n".join(result)


def generate_communication(target: Target, grid: Grid):


    communication_name = "communication"

    result = []

    add_packinfo = "{pi_creation}\ncommunication->addPackInfo( {name} );"

    if target == Target.GPU:
        if grid == Grid.UNIFORM:
            result.append("bool cudaEnabledMPI = parameters.getParameter<bool>( \"cudaEnabledMPI\", false );")
            result.append("if(walberla::mpi::MPIManager::instance()->numProcesses() == 1)")
            result.append("    cudaEnabledMPI = false;\n")

            scheme            = f"auto {communication_name} = std::make_shared<walberla::gpu::communication::UniformGPUScheme<CommunicationStencil_T>>(blocks, cudaEnabledMPI);"
            pdf_packinfo      = f"auto pdfPackInfo = std::make_shared< walberla::lbm_generated::UniformGeneratedGPUPdfPackInfo<{get_field_type(Fields.PDF, target)}> >({get_field_variable(Fields.PDF, target)});"
            packinfo_density  = f"auto densityPackInfo = std::make_shared< walberla::gpu::communication::MemcpyPackInfo<GPUField_T<real_t>> >({get_field_variable(Fields.DENSITY, target)});"
            packinfo_velocity = f"auto velocityPackInfo = std::make_shared< walberla::gpu::communication::MemcpyPackInfo<GPUField_T<real_t>> >({get_field_variable(Fields.VELOCITY, target)});"
        else:
            raise ValueError("Non-uniform GPU is currently not supported")
    elif target == Target.CPU:
        if grid == Grid.UNIFORM:
            scheme            = f"auto {communication_name} = std::make_shared<walberla::blockforest::communication::UniformBufferedScheme<CommunicationStencil_T>>(blocks);"
            pdf_packinfo      = f"auto pdfPackInfo = std::make_shared< walberla::lbm_generated::UniformGeneratedPdfPackInfo<{get_field_type(Fields.PDF)}> >({get_field_variable(Fields.PDF, target)});"
            packinfo_density  = f"auto densityPackInfo = std::make_shared< walberla::field::communication::PackInfo<{get_field_type(Fields.DENSITY)}> >({get_field_variable(Fields.DENSITY, target)});"
            packinfo_velocity = f"auto velocityPackInfo = std::make_shared< walberla::field::communication::PackInfo<{get_field_type(Fields.VELOCITY)}> >({get_field_variable(Fields.VELOCITY, target)});"
        else:
            scheme            = f"auto {communication_name} = std::make_shared<walberla::blockforest::communication::NonUniformBufferedScheme<CommunicationStencil_T>>(blocks);"
            pdf_packinfo      = f"auto pdfPackInfo = walberla::lbm_generated::setupNonuniformPdfCommunication<{get_field_type(Fields.PDF)}>(blocks, {get_field_variable(Fields.PDF, target)});"
            packinfo_density  = f"auto densityPackInfo = std::make_shared< walberla::field::refinement::PackInfo<{get_field_type(Fields.DENSITY)}, CommunicationStencil_T> >({get_field_variable(Fields.DENSITY, target)});"
            packinfo_velocity = f"auto velocityPackInfo = std::make_shared< walberla::field::refinement::PackInfo<{get_field_type(Fields.VELOCITY)}, CommunicationStencil_T> >({get_field_variable(Fields.VELOCITY, target)});"
    else:
        raise ValueError("Unknown target.")

    result.append(scheme)
    result.append(add_packinfo.format(name="pdfPackInfo", pi_creation=pdf_packinfo))
    result.append(add_packinfo.format(name="densityPackInfo", pi_creation=packinfo_density))
    result.append(add_packinfo.format(name="velocityPackInfo", pi_creation=packinfo_velocity))

    return "\n".join(result)


def generate_tke_output(target: Target, field_map: set[Fields], runtime_gpu_indexing: bool):
    if not Fields.STRAIN_RATE in field_map:
        return ""

    result = []
    result.append(f"walberla::pystencils::waLBerlaWind_StrainRateWriter strainRateWriter("
                  f"{get_field_variable(Fields.FORCE, target)}, {get_field_variable(Fields.OMEGA, target)}, "
                  f"{get_field_variable(Fields.PDF, target)}, {get_field_variable(Fields.STRAIN_RATE, target)}"
                  f"{', gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]' if runtime_gpu_indexing else ''});")

    return "\n".join(result)


def generate_flow_driver(target: Target, runtime_gpu_indexing: bool):

    result = [f"walberla::wind::FlowDriverCollection flowDriver(blocks, &timeloop, globalConfig, domainSetup, "
              f"{get_field_variable(Fields.FORCE, target)}, {get_field_variable(Fields.VELOCITY, target)}, fieldGhostLayers{', gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]' if runtime_gpu_indexing else ''});",
              "for(auto & block : *blocks) { flowDriver(&block); }"]

    return "\n".join(result)



def generate_output_before_functions(target: Target, vtk_output: bool, field_map: set[Fields]):
    if not vtk_output:
        return ""

    tke_output = Fields.EDDY_VISCOSITY in field_map or Fields.STRAIN_RATE in field_map

    if target == Target.CPU and not tke_output:
        return ""

    indent = "    "

    result = []
    if tke_output:
        result.append("turbineOutput.addBeforeFunction( [&]() {")
        result.append("for( auto & block : *blocks ) {")
        if Fields.STRAIN_RATE in field_map:
            result.append("    strainRateWriter(&block);")
        result.append("}});")

    if target == Target.GPU:
        result.append("turbineOutput.addBeforeFunction( [&]() {")
        result.append(indent + device_to_host_field_copy(Fields.DENSITY))
        result.append(indent + device_to_host_field_copy(Fields.VELOCITY))
        result.append(indent + device_to_host_field_copy(Fields.FORCE))
        if vtk_output:
            if Fields.MEAN_VELOCITY_OUTPUT in field_map:
                result.append(indent + device_to_host_field_copy(Fields.MEAN_VELOCITY_OUTPUT))
            if Fields.SUM_OF_SQUARES in field_map:
                result.append(indent + device_to_host_field_copy(Fields.SUM_OF_SQUARES))
            if Fields.SUM_OF_CUBES in field_map:
                result.append(indent + device_to_host_field_copy(Fields.SUM_OF_CUBES))

            if Fields.EDDY_VISCOSITY in field_map:
                result.append(indent + device_to_host_field_copy(Fields.EDDY_VISCOSITY))
            if Fields.STRAIN_RATE in field_map:
                result.append(indent + device_to_host_field_copy(Fields.STRAIN_RATE))

        result.append(f"\n{indent}cudaDeviceSynchronize();")

        result.append("});")

    return "\n".join(result)


def generate_stability_checker():
    return (f"walberla::makeSharedFunctor( walberla::field::makeStabilityChecker< {get_field_type(Fields.PDF)}, "
            f"{get_field_type(Fields.FLAG)} >(globalConfig, blocks, {get_field_variable(Fields.PDF, Target.CPU)}, {get_field_variable(Fields.FLAG)}, FluidFlagUID ))")


def generate_performance_evaluation():
    result = [
        f"walberla::lbm::PerformanceEvaluation<{get_field_type(Fields.FLAG)}> performance(blocks, {get_field_variable(Fields.FLAG)}, FluidFlagUID);",
        "performance.logResultOnRoot(timesteps, time);"
    ]

    return "\n".join(result)


def generate_field_output_list(field_map: set[Fields]):
    result = []

    # order set to have consistent output
    field_list = list(field_map)
    field_list = sorted(field_list, key=lambda field: field.name)

    exclude = [Fields.PDF, Fields.FLAG, Fields.MEAN_VELOCITY_WFB]

    for field in field_list:

        # manually exclude some fields that you typically don't want
        if field in exclude:
            continue

        field_name = get_field_name(field).replace("field", "").title().replace(" ", "") + ';'
        result.append(field_name)

    return "\n".join(result)


def set_additional_headers(target: Target, grid: Grid):
    walberla_headers = ["field/Layout.h", "lbm_generated/field/PdfField.h", "lbm_generated/field/AddToStorage.h"]
    turbine_headers = []

    if target == Target.GPU:
        walberla_headers.append("lbm_generated/gpu/GPUPdfField.h")
        walberla_headers.append("lbm_generated/gpu/AddToStorage.h")
        walberla_headers.append("gpu/AddGPUFieldToStorage.h")
        walberla_headers.append("gpu/HostFieldAllocator.h")
        walberla_headers.append("gpu/communication/MemcpyPackInfo.h")
        walberla_headers.append("gpu/ShiftedPeriodicity.h")
        turbine_headers.append("turbine_topology/gpu/Tree.h")

        if grid == Grid.UNIFORM:
            walberla_headers.append("lbm_generated/gpu/UniformGeneratedGPUPdfPackInfo.h")
            walberla_headers.append("gpu/communication/UniformGPUScheme.h")
        else:
            walberla_headers.append("lbm_generated/gpu/NonuniformGeneratedGPUPdfPackInfo.h")
            walberla_headers.append("gpu/communication/NonUniformGPUScheme.h")

    elif target == Target.CPU:
        walberla_headers.append("field/allocation/FieldAllocator.h")
        walberla_headers.append("boundary/ShiftedPeriodicity.h")
        turbine_headers.append("turbine_topology/cpu/Tree.h")

        if grid == Grid.UNIFORM:
            walberla_headers.append("blockforest/communication/UniformBufferedScheme.h")
            walberla_headers.append("lbm_generated/communication/UniformGeneratedPdfPackInfo.h")
        else:
            walberla_headers.append("blockforest/communication/NonUniformBufferedScheme.h")
            walberla_headers.append("lbm_generated/communication/NonuniformGeneratedPdfPackInfo.h")
    else:
        raise ValueError("Invalid target")

    if grid != Grid.UNIFORM:
        turbine_headers.append("walberla_helper/refinement/CustomRecursiveTimeStep.h")

    result = ["additional_headers = {"]

    headers = walberla_headers + turbine_headers

    for header in headers:
        result.append(f"\t\"{header}\",")

    result.append("}")

    return "\n".join(result)


def generate_shifted_periodicity(target: Target):
    result = ["const int periodicShiftValue = boundariesConfig.getParameter<int>(\"periodicShiftValue\", int(real_t(0.33) * real_t(domainSetup.domainSize_[2])));"]

    namespace = "boundary" if target == Target.CPU else "gpu"
    pdf_field_type = get_field_type(Fields.PDF) if target == Target.CPU else get_field_type(Fields.PDF, target)

    result.append(f"walberla::{namespace}::ShiftedPeriodicity{'GPU' if target == Target.GPU else ''}<{pdf_field_type}> shiftedPeriodicity (blocks, {get_field_variable(Fields.PDF, target)}, fieldGhostLayers, 0, 1, periodicShiftValue);")
    return "\n".join(result)


def add_turbine_filters_to_jinja_env(jinja_env, jinja_globals):
    jinja_globals['Stencil'] = Stencil
    jinja_globals['Target'] = Target
    jinja_globals['Grid'] = Grid
    jinja_globals['Fields'] = Fields
    jinja_globals['FieldType'] = FieldType
    jinja_globals['StreamingPattern'] = StreamingPattern
    jinja_env.filters['generate_gpu_field_creation'] = generate_gpu_field_creation
    jinja_env.filters['generate_field_creation'] = generate_field_creation
    jinja_env.filters['generate_field_initialisation'] = generate_field_initialisation
    jinja_env.filters['generate_snapshot_creation'] = generate_snapshot_creation
    jinja_env.filters['generate_sweep_collection'] = generate_sweep_collection
    jinja_env.filters['generate_boundary_collection'] = generate_boundary_collection
    jinja_env.filters['generate_set_field_ids'] = generate_set_field_ids
    jinja_env.filters['generate_field_resetter'] = generate_field_resetter
    jinja_env.filters['generate_communication'] = generate_communication
    jinja_env.filters['generate_tke_output'] = generate_tke_output
    jinja_env.filters['generate_output_before_functions'] = generate_output_before_functions
    jinja_env.filters['generate_shifted_periodicity'] = generate_shifted_periodicity
    jinja_env.filters['generate_flow_driver'] = generate_flow_driver
    jinja_globals['generate_stability_checker'] = generate_stability_checker
    jinja_globals['Welford'] = Welford
    jinja_globals['Timeloop'] = Timeloop
    jinja_globals['generate_performance_evaluation'] = generate_performance_evaluation
    jinja_env.filters['generate_field_output_list'] = generate_field_output_list
    jinja_env.filters['set_additional_headers'] = set_additional_headers
