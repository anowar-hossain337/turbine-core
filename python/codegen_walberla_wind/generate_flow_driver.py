import pystencils as ps
import lbmpy
import numpy as np

from jinja2 import Environment, PackageLoader, StrictUndefined

from pystencils_walberla.cmake_integration import CodeGenerationContext
from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
from pystencils_walberla import generate_sweep, generate_sweep_collection, function_generator

from .flow_drivers import *


def generate_flow_driver_collection(generation_context: CodeGenerationContext, class_name: str, namespace: str,
                                    target: ps.Target, force_field: [ps.Field, ps.Field.Access],
                                    velocity_field: [ps.Field, ps.Field.Access],
                                    ghost_layers_to_include=0, gpu_indexing_params=None, max_threads=None,
                                    **create_kernel_parameters):

    #TODO add option not to create ALL drivers but only a selection

    includes = []
    driver_class_dict = {}
    object_names = []
    kernels = []

    driver_kernels = []

    # create flow driver classes
    for driver_enum, driver in set(driver_dict.items()):

        driver_class_name = driver.__name__ + "Driver"
        driver_obj = driver(force_field=force_field, velocity_field=velocity_field)

        functor = function_generator(ctx=generation_context, class_name=driver_class_name, namespace=namespace,
                                     assignments=driver_obj(), ghost_layers_to_include=ghost_layers_to_include,
                                     target=target, gpu_indexing_params=gpu_indexing_params, max_threads=max_threads,
                                     **create_kernel_parameters)

        context = functor()

        driver_kernels.append(functor)

        includes.append(f"\"{driver_class_name}.h\"")
        driver_class_dict[driver_enum] = context['function_name']
        object_names.append(f"{context['function_name'][0].lower() + context['function_name'][1:]}Object")
        kernels.append(context['kernel'])

    collection_name = "Drivers"
    generate_sweep_collection(generation_context, class_name=collection_name, function_generators=driver_kernels)

    parameter_list = []
    parameters_to_skip = []
    # First field pointer
    for kernel_info in kernels:
        for param in kernel_info.parameters:
            if param.is_field_pointer and param.field_name not in parameters_to_skip:
                parameter_list.append(param.field_name + "FieldID")
                parameters_to_skip.append(param.field_name)

    # Then free parameters
    for kernel_info in kernels:
        for param in kernel_info.parameters:
            if not param.is_field_parameter and param.symbol.name not in parameters_to_skip:
                parameter_list.append(param.symbol.name)
                parameters_to_skip.append(param.symbol.name)

    # create flow driver collection

    jinja_context = {
        'kernel_list': kernels,
        'class_name': class_name,
        'collection_name': collection_name,
        'collection_object_name': collection_name[0].lower() + collection_name[1:] + "Object_",
        'parameter_list': ", ".join(parameter_list),
        'target': target,
        'namespace': namespace,
        'includes': includes,
        'driver_classes': driver_class_dict,
        'object_names': object_names,
        'data_type': 'double' if generation_context.double_accuracy else 'float',
        'runtime_gpu_indexing': np.any(['gpuBlockSize' in param for param in parameter_list]) #FIXME: only works with our definition of indices, not in the general case
    }

    env = Environment(loader=PackageLoader('codegen_walberla_wind'), undefined=StrictUndefined)
    env.globals.update(zip=zip)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template("FlowDriverCollection.tmpl.h").render(**jinja_context)

    generation_context.write_file(f"{class_name}.h", header)
