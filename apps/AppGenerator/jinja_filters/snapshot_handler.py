from .enums import *
from .welford import Welford


def generate_snapshot_creation(generate_wfb: bool, vtk_output: bool,
                               target: Target, field_map: set[Fields], indent: int,
                               welford_wfb: Welford, welford_output: Welford,
                               timeloop_name: str):

    indent = " " * indent
    indent += "      "

    create_snapshot_string = (indent + "   " + "domainSetup.snapshotHandler_->createSnapshot(timeloop.getCurrentTimeStep(), blocks, {pdf_field}, {force_field}"
                                               "{optional_arguments});")
    result = [" [&](){"]
    result.append(indent + f"if ((domainSetup.snapshotHandler_->storeSnapshot) && "
                            f"({timeloop_name}.getCurrentTimeStep() % domainSetup.snapshotHandler_->snapshotFrequency == 0) && "
                            f"({timeloop_name}.getCurrentTimeStep() > 0)) {{")

    if target == Target.GPU:
        # field copies
        result.append(indent + "   " + device_to_host_field_copy(Fields.PDF))
        result.append(indent + "   " + device_to_host_field_copy(Fields.FORCE))

        if generate_wfb and Fields.MEAN_VELOCITY_WFB in field_map:
            result.append(indent + "   " + device_to_host_field_copy(Fields.MEAN_VELOCITY_WFB))

        if vtk_output:
            if not generate_wfb and Fields.MEAN_VELOCITY_OUTPUT in field_map:
                result.append(indent + "   " + device_to_host_field_copy(Fields.MEAN_VELOCITY_OUTPUT))
            if Fields.SUM_OF_SQUARES in field_map:
                result.append(indent + "   " + device_to_host_field_copy(Fields.SUM_OF_SQUARES))
            if Fields.SUM_OF_CUBES in field_map:
                result.append(indent + "   " + device_to_host_field_copy(Fields.SUM_OF_CUBES))

        result.append(indent + "   gpuDeviceSynchronize();")

    indent_string = indent + " " * len("   domainSetup.snapshotHandler_->createSnapshot(")

    optional_arguments = []
    if generate_wfb or vtk_output:

        assert (welford_wfb.mean_field in field_map) or (welford_output.mean_field in field_map), \
            "At least one mean velocity field must be provided in the field map."

        if welford_wfb.active:
            optional_arguments.append(f",\n{indent_string}uint_t({welford_wfb.sweep_name()}.getCounter())")
            optional_arguments.append(f"{get_field_variable(welford_wfb.mean_field, Target.CPU)}")
        else:
            optional_arguments.append(f",\n{indent_string}uint_t({welford_output.sweep_name()}.getCounter())")
            optional_arguments.append(f"{get_field_variable(welford_output.mean_field, Target.CPU)}")


        if len(welford_output.welford_fields) >= 3:
            optional_arguments.append(f"{get_field_variable(welford_output.sos_field, Target.CPU)}")
        if len(welford_output.welford_fields) >= 4:
            optional_arguments.append(f"{get_field_variable(welford_output.soc_field, Target.CPU)}")

    result.append(create_snapshot_string.format(pdf_field=get_field_variable(Fields.PDF, Target.CPU),
                                                force_field=get_field_variable(Fields.FORCE, Target.CPU),
                                                optional_arguments=", ".join(optional_arguments)))

    result.append(indent + "}")
    result.append(indent + "}")

    return "\n".join(result)