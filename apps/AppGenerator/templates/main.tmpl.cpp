{% set welford_wfb = Welford(welford_wfb_name, target, generate_wfb,
                             field=Fields.VELOCITY, mean_field=Fields.MEAN_VELOCITY_WFB) -%}
{% set welford_output = Welford(welford_output_name, target, vtk_output and Fields.MEAN_VELOCITY_OUTPUT in field_map,
                                field=Fields.VELOCITY, mean_field=Fields.MEAN_VELOCITY_OUTPUT,
                                sos_field=Fields.SUM_OF_SQUARES if Fields.SUM_OF_SQUARES in field_map else None,
                                soc_field=Fields.SUM_OF_CUBES if Fields.SUM_OF_CUBES in field_map else None) -%}
{% set welford_eddy_viscosity = Welford(welford_eddy_viscosity_name, target, vtk_output and Fields.MEAN_EDDY_VISCOSITY in field_map,
                                        field=Fields.EDDY_VISCOSITY, mean_field=Fields.MEAN_EDDY_VISCOSITY) -%}
{% set welford_strain_rate = Welford(welford_strain_rate_name, target, vtk_output and Fields.MEAN_STRAIN_RATE in field_map,
                                     field=Fields.STRAIN_RATE, mean_field=Fields.MEAN_STRAIN_RATE) -%}

{% set timeloop = Timeloop("timeloop", grid, target, field_map, vtk_output, with_snapshot, generate_wfb,
                           turbine_communication, generate_control, welford_wfb, welford_output,
                           welford_eddy_viscosity, welford_strain_rate) -%}

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>

#include <blockforest/Initialization.h>
#include <core/all.h>
#include <core/Environment.h>
#include <core/mpi/all.h>
#include <domain_decomposition/all.h>
#include <field/all.h>
#include <geometry/all.h>
#include <lbm/all.h>
#include <lbm/vtk/all.h>
#include <lbm_generated/field/AddToStorage.h>
#include <timeloop/all.h>

#include "conversion/Conversion.h"
#include "domain/BoundaryHandling.h"
#include "domain/BoundarySetter.h"
#include "domain/BoundarySetup.h"
#include "domain/DomainInitialiser.h"
#include "domain/DomainSetup.h"
#include "farm_creator/WindFarm.h"
#include "farm_creator/ConfigTurbineCreator.h"
#include "output/Logging.h"
#include "output/TurbineOutput.h"
#include "walberla_helper/field/all.h"
#include "wind_turbine_core/ProjectDefines.h"

#include "waLBerlaWind_KernelInfo.h"

namespace turbine_core {

    static const uint_t fieldGhostLayers = {% if grid == Grid.UNIFORM %} 1 {% else %} 2 {% endif -%};

    // define different interpolators and distributors

    template<class Kernel_T>
    using KernelFieldInterpolator_T = projectors::KernelFieldInterpolator<field::Field<real_t>, Kernel_T>;

    template<class Kernel_T>
    using KernelDistributor_T = projectors::KernelFieldDistributor<field::Field<real_t>, Kernel_T>;

    using InterpolationKernel_T = projectors::SmoothedDiracDeltaKernel;
    using DistributionKernel_T = projectors::SmoothedDiracDeltaKernel;

    using WindFarm_T = WindFarm<topology::{{target.name.lower()}}::Tree, creator::ConfigTurbineCreator,
            ScalarField_T, VectorField_T, VectorField_T,
            KernelFieldInterpolator_T<InterpolationKernel_T>,
            KernelDistributor_T<DistributionKernel_T>>;

    using TurbineOutput_T = output::TurbineOutput<
            FlagField_T, ScalarField_T, VectorField_T, SecondOrderTensorField_T, ThirdOrderTensorField_T,
            PdfField_T, Stencil_T, StorageSpecification_T::zeroCenteredPDFs, StorageSpecification_T::compressible
    >;

    {% if target == Target.GPU -%}
    template<typename Type_T>
    using GPUField_T = walberla::gpu::GPUField<Type_T>;
    {%- endif %}

    int main(int argc, char** argv) {

        walberla::Environment walberlaEnv(argc, argv);

        walberla::logging::Logging::instance()->setLogLevel(walberla::logging::Logging::INFO);

        /// general simulation parameters
        auto globalConfig = walberlaEnv.config();
        auto parameters = globalConfig->getOneBlock("Parameters");

        uint_t timesteps = parameters.getParameter<uint_t>("timesteps", uint_t(10)); ++timesteps;
        const real_t remainingTimeLoggerFrequency = parameters.getParameter<real_t>("remainingTimeLoggerFrequency", real_t(3.0)); // in seconds

        const real_t length_SI = parameters.getParameter<real_t>("length_SI");
        const real_t length_LU = parameters.getParameter<real_t>("length_LU");
        const real_t velocity_SI = parameters.getParameter<real_t>("velocity_SI");
        const real_t velocity_LU = parameters.getParameter<real_t>("velocity_LU");
        const real_t viscosity_SI = parameters.getParameter<real_t>("viscosity_SI");
        const real_t density_SI = parameters.getParameter<real_t>("density_SI");

        Conversion::calculateConversionFactors(length_SI, length_LU, velocity_SI, velocity_LU, density_SI);
        Conversion::print();

//        projectors::GaussianFunction5<>::setSigma(real_t(7.07106781));

        /// parameter conversion
        real_t viscosity_LU = (viscosity_SI / density_SI) * Conversion::C_t() / Conversion::C_l() / Conversion::C_l();

        WALBERLA_CHECK_GREATER(viscosity_LU, real_t(0.0), "Negative LU velocity - Did you calculate the conversion factors?")

        // calculate omega on coarsest level
        const real_t omega = walberla::lbm::collision_model::omegaFromViscosity(viscosity_LU);
        WALBERLA_LOG_INFO_ON_ROOT("omega = " << omega)

        /// TURBINE

        auto windFarmConfig = globalConfig->getBlock("WindFarm");
        {% if turbine_output -%}
        const std::string turbineOutputDirectory = windFarmConfig.getParameter<std::string>("outputFolder", "");
        const uint_t forceOutputStart = windFarmConfig.getParameter<uint_t>("startingTimeStepForces", 0);
        const uint_t forceOutputFrequency = windFarmConfig.getParameter<uint_t>("frequencyOutputForces", 1);
        {%- else %}
        const std::string turbineOutputDirectory{};
        {%- endif %}

        WindFarm_T farm(turbineOutputDirectory, windFarmConfig);
        farm.createGeometry();
        farm.calculateTurbineAABBs();

        /// BLOCKFOREST

        WALBERLA_LOG_INFO_ON_ROOT("Creating blockforest...")

        auto boundariesConfig = globalConfig->getOneBlock( "Boundaries" );
        domain::BoundarySetup boundarySetup(boundariesConfig);

        domain::DomainSetup domainSetup ( globalConfig, boundarySetup.periodicity() );
        {% if grid == Grid.UNIFORM -%}
        auto blocks = domainSetup.createUniformBlockForest(farm.getTurbineAABBs());
        {%- else %}
        WALBERLA_ASSERT(!domainSetup.snapshotHandler_->loadSnapshot, "Restart is currently not supported for non-uniform grids.")

        refinement::StaticRefinementHandler<WindFarm_T> staticRefinementHandler{ globalConfig->getOneBlock("DomainSetup"), &farm, boundarySetup.environmentSetup() };
        auto blocks = domainSetup.createStructuredBlockForest(staticRefinementHandler, fieldGhostLayers);
        {%- endif %}

        farm.setBlockForest(blocks);
        farm.initialiseForceAndControlModel();

        /// FIELDS

        WALBERLA_LOG_INFO_ON_ROOT("Creating fields...")

        //layout - DO NOT CHANGE
        auto layout = codegen::KernelInfo::layout;
        const StorageSpecification_T storageSpec{};

        {{field_map|generate_field_creation(target,with_snapshot,generate_wfb,vtk_output)|indent(8)}}

        /// ABL calculations

        auto initParams = globalConfig->getOneBlock("Initialisation");
        const auto initType = domain::DomainInitialisation::toType(initParams.getParameter<std::string>("type"));
        const Vector3<real_t> initialVelocity = initParams.getParameter<Vector3<real_t>>("initialVelocity");

        const real_t kappa = parameters.getParameter<real_t>("kappa", real_t(0.42));
        const real_t B = parameters.getParameter<real_t>("B", real_t(5.5));
        const real_t roughnessLengthRatio = parameters.getParameter<real_t>("roughnessLengthRatio", real_t(1e-4));
        const real_t referenceHeight = parameters.getParameter<real_t>("referenceHeight_LU", real_t(-1));
        {% if generate_wfb -%}
        const uint32_t samplingHeight = parameters.getParameter<uint32_t>("samplingHeight_LU", uint32_t(0));
        {% endif -%}
        const real_t roughnessLength = roughnessLengthRatio * referenceHeight;

        {% if vtk_output -%}
        {{Welford.read_interval()|indent(8)}}
        {% endif -%}

        // uTau from MOST
        const real_t uTau = kappa * initialVelocity[0] / std::log(real_t(1) / roughnessLengthRatio);

        /// DOMAIN INITIALISATION

        WALBERLA_LOG_INFO_ON_ROOT("Initialise domain...")
        {{field_map|generate_field_initialisation(with_snapshot, vtk_output, generate_wfb, streaming_pattern)|indent(8)}}

        {% if target == Target.GPU -%}
        {{field_map|generate_gpu_field_creation|indent(8)}}
        {% endif -%}

        {% if target == Target.GPU and runtime_gpu_indexing -%}
        Vector3 <int32_t> gpuBlockSize = parameters.getParameter<Vector3<int32_t>>("gpuBlockSize", Vector3<int32_t>(256, 1, 1));
        {% endif -%}
        {{target|generate_sweep_collection(field_map, runtime_gpu_indexing)|indent(8)}}

        /// BOUNDARIES

        WALBERLA_LOG_INFO_ON_ROOT("Setting up boundaries...")

        boundarySetup.fillFlagFieldFromConfig<FlagField_T>(blocks, flagFieldID, FluidFlagUID,
                                                           NoSlipFlagUID, WFBFlagUID, SymmetryFlagUID, UniformInflowFlagUID,
                                                           LogLawInflowFlagUID, OutflowFlagUID);

        {% if generate_wfb -%}
        if(boundarySetup.wallType() == WallSetup::WFB){
            WALBERLA_CHECK(referenceHeight > real_t(0), "Reference height must be given in the parameter file for wall function boundary conditions.")
            WALBERLA_CHECK(boundarySetup.environmentSetup() != EnvironmentSetup::Tunnel,
                           "Wall-function bounce is currently not supported for the tunnel environment")
            WALBERLA_CHECK(samplingHeight > uint_t(0), "Sampling height must be given in the parameter file for wall function boundary conditions.")
        }
        {% endif %}

        {{target|generate_boundary_collection(generate_wfb)|indent(8)}}

        {{target|generate_shifted_periodicity|indent(8)}}

        /// set IDs for farm
        {{target|generate_set_field_ids|indent(8)}}

        /// field resetter
        {{target|generate_field_resetter(welford_output, runtime_gpu_indexing)|indent(8)}}

        {{welford_wfb.object_creation(runtime_gpu_indexing)|indent(8)}}
        {{welford_output.object_creation(runtime_gpu_indexing)|indent(8)}}
        {{welford_eddy_viscosity.object_creation(runtime_gpu_indexing)|indent(8)}}
        {{welford_strain_rate.object_creation(runtime_gpu_indexing)|indent(8)}}

        /// COMMUNICATION

        WALBERLA_LOG_INFO_ON_ROOT("Set up communication...")

        // create communication for fields
        {{target|generate_communication(grid)|indent(8)}}

        /// TIMELOOP

        WALBERLA_LOG_INFO_ON_ROOT("Creating time loop...")

        {% if target == Target.CPU -%}
        using TimingPolicy_T = walberla::timing::WcPolicy;
        {% elif target == Target.GPU -%}
#ifdef NDEBUG
        using TimingPolicy_T = walberla::timing::WcPolicy;
#else
        using TimingPolicy_T = walberla::timing::DeviceSynchronizePolicy;
#endif
        {% endif -%}

        // create time loop
        {{timeloop.declaration()|indent(8)}}

        farm.setTimeloop(&{{timeloop.name}});

        /// driving forces
        {{target|generate_flow_driver(runtime_gpu_indexing)|indent(8)}}

        uint_t maxKernelWidth = 1;
        if constexpr (std::is_same_v<WindFarm_T::Interpolator, KernelFieldInterpolator_T<InterpolationKernel_T>>) {
            maxKernelWidth = std::max({maxKernelWidth, InterpolationKernel_T::xWidth, InterpolationKernel_T::yWidth, InterpolationKernel_T::zWidth});
        }
        if constexpr (std::is_same_v<WindFarm_T::Distributor, KernelDistributor_T<DistributionKernel_T >>) {
            maxKernelWidth = std::max({maxKernelWidth, DistributionKernel_T::xWidth, DistributionKernel_T::yWidth, DistributionKernel_T::zWidth});
        }

        walberla::Config::Blocks turbineBlocks{};
        windFarmConfig.getBlocks(turbineBlocks);
        auto numTurbines = turbineBlocks.size();
        uint_t diameter = 0;
        //TODO will we ever have different models during the simulation ?
        for (auto & turbineConfig : turbineBlocks) {
            const std::string aeroModel = turbineConfig.getParameter<std::string>( "aeroModel" );
            if(aeroModel == "rotorDisk") {
                diameter = math::max(diameter, turbineConfig.getParameter<uint_t>("diameter_LU"));
            }
        }

        auto evaluateDensityAndVelocity = [&farm](walberla::IBlock * block){farm.evaluateDensityAndVelocity(block);};
        {% if generate_control -%}
        auto applyControl = [&farm](){farm.applyControl();};
        {% endif -%}
        {% if turbine_communication -%}
        auto syncForces = [&](){farm.syncNextNeighbour(mpi::SYNC_FORCE, maxKernelWidth+uint_t(real_t(diameter)/real_t(2.0)), false);};
        auto syncMacros = [&](){farm.syncNextNeighbour(mpi::SYNC_MACRO, 1, false);};
        {% endif %}
        auto spreadForces = [&farm](walberla::IBlock * block){farm.spreadForces(block);};

        // add LBM sweep and communication to time loop
        // output
        TurbineOutput_T turbineOutput{
                blocks, &{{timeloop.name}},
                globalConfig->getOneBlock("Output"), FluidFlagUID,
                fieldMap
        };

        {{target|generate_tke_output(field_map, runtime_gpu_indexing)|indent(8)}}

        {{target|generate_output_before_functions(vtk_output,field_map)|indent(8)}}

        {% if target == Target.CPU -%}
        auto output = [&turbineOutput]() { turbineOutput.write(); };
        {% elif target == Target.GPU -%}
        auto output = [&turbineOutput](){ cudaDeviceSynchronize(); turbineOutput.write(); };
        {% endif %}

        {% if welford_wfb.active -%}
        {{welford_wfb.create_lambda()|indent(8)}}
        {%- endif %}
        {% if welford_output.active -%}
        {{welford_output.create_lambda()|indent(8)}}
        {%- endif %}
        {% if welford_eddy_viscosity.active -%}
        {{welford_eddy_viscosity.create_lambda()|indent(8)}}
        {%- endif %}
        {% if welford_strain_rate.active -%}
        {{welford_strain_rate.create_lambda()|indent(8)}}
        {%- endif %}

        {{timeloop.add_to_timeloop()|indent(8)}}

        // LBM stability check
        {{timeloop.add_after_timestep(function=("LBM stability check", generate_stability_checker()))|indent(8)}}

        // log remaining time
        {{timeloop.add_after_timestep(function=("Remaining time logger", timeloop.remaining_time_logger()))|indent(8)}}

        WALBERLA_LOG_INFO_ON_ROOT("Running timeloop...")

        WALBERLA_MPI_WORLD_BARRIER()

        walberla::timing::TimingPool<TimingPolicy_T> timing;
        {% if turbine_output -%}
        timing.registerTimer("Force Output");
        {% endif %}

        walberla::WcTimer timer;
        timer.start();

        for( uint_t i = 0; i < timesteps; ++i ) {

            {{timeloop.single_step()}}

            {% if turbine_output -%}
            timing["Force Output"].start();
            // farm.callback(Component::Output::ORIENTATIONS);
            for(auto & block : *blocks) {
                farm.writeForceOutput(&block, forceOutputStart, forceOutputFrequency);
            }
            farm.writeTurbinePowerAndThrust();
            timing["Force Output"].end();
            {% endif %}
        }

        timer.end();

        double time = timer.max();
        walberla::mpi::reduceInplace(time, walberla::mpi::MAX);

        const auto timeloopTiming = timing.getReduced();
        WALBERLA_LOG_INFO_ON_ROOT("Timeloop timing:\n" << *timeloopTiming)

        {{generate_performance_evaluation()|indent(8)}}

        farm.calculatePowerAndThrust();

        logging::WindTurbineApplicationLogger<domain::DomainSetup, domain::BoundarySetup, codegen::KernelInfo, TimingPolicy_T> logger {
            domainSetup, boundarySetup, windFarmConfig, &timing, performance.loggingString(timesteps, time)
        };

        return EXIT_SUCCESS;

    } // main

} // namespace turbine_core

int main(int argc, char** argv) {
    return turbine_core::main(argc, argv);
}
