#include <iostream>
#include <cmath>
#include <vector>
#include <functional>

#include <core/Environment.h>

#include "output/Logging.h"

//TODO
#include <blockforest/Initialization.h>
#include <boundary/ShiftedPeriodicity.h>
#include <core/all.h>
#include <core/mpi/all.h>
#include <domain_decomposition/all.h>
#include <field/all.h>
#include <geometry/all.h>
#include <lbm_generated/field/AddToStorage.h>
#include <lbm_generated/communication/UniformGeneratedPdfPackInfo.h>
#include <lbm/all.h>
#include <lbm/vtk/all.h>
#include <timeloop/all.h>

#include "walberla_helper/field/all.h"

#include "domain/BoundarySetup.h"
#include "domain/BoundarySetter.h"
#include "domain/DomainSetup.h"
#include "domain/DomainInitialiser.h"

#include "farm_creator/WindFarm.h"
#include "farm_creator/ConfigTurbineCreator.h"

#include "wind_turbine_core/ProjectDefines.h"

#include "conversion/Conversion.h"

#include "turbine_topology/cpu/Tree.h"

#include "WindTurbineDefines.h"

#include "output/TurbineOutput.h"

#include "waLBerlaWind_EddyViscosityWriter.h"
#include "waLBerlaWind_StrainRateWriter.h"
#include "waLBerlaWind_Welford.h"
#include "waLBerlaWind_WelfordSOP.h"
#include "waLBerlaWind_MacroGetter.h"
#include "domain/BoundaryHandling.h"

namespace turbine_core {

    using InterpolationKernel_T = projectors::SmoothedDiracDeltaKernel;
    using DistributionKernel_T = projectors::SmoothedDiracDeltaKernel;

    using WindFarm_T = WindFarm<topology::cpu::Tree, creator::ConfigTurbineCreator,
            ScalarField_T, VectorField_T, VectorField_T,
//            TrilinearMacroscopicFieldInterpolator,
            KernelFieldInterpolator_T<InterpolationKernel_T>,
            KernelDistributor_T<DistributionKernel_T>>;

    using TurbineOutput_T = output::TurbineOutput<
            FlagField_T, ScalarField_T, VectorField_T, SecondOrderTensorField_T, ThirdOrderTensorField_T,
            PdfField_T, Stencil_T, StorageSpecification_T::zeroCenteredPDFs, StorageSpecification_T::compressible
    >;

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

        projectors::GaussianFunction5<>::setSigma(real_t(7.07106781));

        /// parameter conversion
        real_t viscosity_LU = (viscosity_SI / density_SI) * Conversion::C_t() / Conversion::C_l() / Conversion::C_l();

        WALBERLA_CHECK_GREATER(viscosity_LU, real_t(0.0), "Negative LU velocity - Did you calculate the conversion factors?")

        // calculate omega on coarsest level
        const real_t omega = walberla::lbm::collision_model::omegaFromViscosity(viscosity_LU);
        WALBERLA_LOG_INFO_ON_ROOT("omega = " << omega)

        /// TURBINE

        auto windFarmConfig = globalConfig->getBlock("WindFarm");
        const std::string turbineOutputDirectory = windFarmConfig.getParameter<std::string>("outputFolder", "");
        const uint_t forceOutputStart = windFarmConfig.getParameter<uint_t>("startingTimeStepForces", 0);
        const uint_t forceOutputFrequency = windFarmConfig.getParameter<uint_t>("frequencyOutputForces", 1);

        WindFarm_T farm(turbineOutputDirectory, windFarmConfig);
        farm.createGeometry();
        farm.calculateTurbineAABBs();

        /// BLOCKFOREST

        WALBERLA_LOG_INFO_ON_ROOT("Creating blockforest...")

        auto boundariesConfig = globalConfig->getOneBlock( "Boundaries" );
        domain::BoundarySetup boundarySetup(boundariesConfig);

        domain::DomainSetup domainSetup ( globalConfig, boundarySetup.periodicity() );
        auto blocks = domainSetup.createUniformBlockForest(farm.getTurbineAABBs());
        farm.setBlockForest(blocks);

        farm.initialiseForceAndControlModel();

        /// FIELDS

        WALBERLA_LOG_INFO_ON_ROOT("Creating fields...")

        std::map<output::Fields::Types, BlockDataID> fieldMap {};

        //layout - DO NOT CHANGE
        auto layout = codegen::KernelInfo::layout;

        BlockDataID pdfFieldID;
        BlockDataID forceFieldID;
        BlockDataID meanVelocityFieldWFBID;
        BlockDataID meanVelocityFieldOutputID;
        BlockDataID sumOfSquaresFieldID;

        BlockDataID nutFieldID = walberla::field::addToStorage<ScalarField_T>(blocks, "eddy viscosity field", real_t(0), layout, fieldGhostLayers);
        BlockDataID strainRateFieldID = walberla::field::addToStorage<SecondOrderTensorField_T>(blocks, "strain rate field", real_t(0), layout, fieldGhostLayers);
        BlockDataID velocityFieldID = walberla::field::addToStorage<VectorField_T>(blocks, "velocity field", real_t(0), layout, fieldGhostLayers);
        BlockDataID densityFieldID = walberla::field::addToStorage<ScalarField_T>(blocks, "density field", real_t(1), layout, fieldGhostLayers);

        BlockDataID flagFieldID = walberla::field::addFlagFieldToStorage<FlagField_T>(blocks, "flag field", fieldGhostLayers);

        const StorageSpecification_T storageSpec{};
        auto allocator = std::make_shared< walberla::field::AllocateAligned<real_t, 64> >();

        uint_t initialWelfordCounter = 0;
        if(domainSetup.snapshotHandler_->loadSnapshot) {
            auto & ssh = domainSetup.snapshotHandler_;

            auto vectorFieldDataHandling = std::make_shared<field::FlattenedFieldHandling<VectorField_T>> (blocks, fieldGhostLayers, layout, allocator);
            const auto forceFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->forceFile).c_str());

            forceFieldID = blocks->loadBlockData(forceFileName, vectorFieldDataHandling, "force field");

            auto pdfFieldDataHandling = std::make_shared< walberla::lbm_generated::internal::PdfFieldHandling< StorageSpecification_T > >( blocks, storageSpec, fieldGhostLayers, layout );

            const auto pdfFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->pdfFile).c_str());
            pdfFieldID = blocks->loadBlockData(pdfFileName, pdfFieldDataHandling, "pdf field");


            if(!ssh->meanVelFile.empty()) {
                initialWelfordCounter = ssh->welfordCounter;
                const auto meanVelFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->meanVelFile).c_str());

                meanVelocityFieldWFBID    = blocks->loadBlockData(meanVelFileName, vectorFieldDataHandling, "mean velocity field WFB");
                meanVelocityFieldOutputID = blocks->loadBlockData(meanVelFileName, vectorFieldDataHandling, "mean velocity field output");

                if(!ssh->sosFile.empty()) {
                    auto secondOrderTensorFieldDataHandling = std::make_shared<field::FlattenedFieldHandling<SecondOrderTensorField_T>> (blocks, fieldGhostLayers, layout, allocator);

                    const auto sosFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->sosFile).c_str());
                    sumOfSquaresFieldID = blocks->loadBlockData(sosFileName, secondOrderTensorFieldDataHandling, "sum of squares field");
                }
            }

        } else {
            forceFieldID = walberla::field::addToStorage<VectorField_T>(blocks, "force field", real_t(0), layout, fieldGhostLayers);

            pdfFieldID = walberla::lbm_generated::addPdfFieldToStorage(blocks, "pdf field", storageSpec, fieldGhostLayers, layout);
            meanVelocityFieldWFBID = walberla::field::addToStorage<VectorField_T>(blocks, "mean velocity field WFB", real_t(0), layout, fieldGhostLayers);
            meanVelocityFieldOutputID = walberla::field::addToStorage<VectorField_T>(blocks, "mean velocity field output", real_t(0), layout, fieldGhostLayers);
            sumOfSquaresFieldID = walberla::field::addToStorage<SecondOrderTensorField_T>(blocks, "sum of squares field", real_t(0), layout, fieldGhostLayers);

        }

        fieldMap[output::Fields::PDF] = pdfFieldID;
        fieldMap[output::Fields::FORCE] = forceFieldID;
        fieldMap[output::Fields::DENSITY] = densityFieldID;
        fieldMap[output::Fields::VELOCITY] = velocityFieldID;
        fieldMap[output::Fields::MEAN_VELOCITY_OUTPUT] = meanVelocityFieldOutputID;
        fieldMap[output::Fields::STRAIN_RATE] = strainRateFieldID;
        fieldMap[output::Fields::EDDY_VISCOSITY] = nutFieldID;
        fieldMap[output::Fields::SUM_OF_SQUARES] = sumOfSquaresFieldID;

        fieldMap[output::Fields::FLAG] = flagFieldID;

        /// ABL calculations

        auto initParams = globalConfig->getOneBlock("Initialisation");
        const auto initType = domain::DomainInitialisation::toType(initParams.getParameter<std::string>("type", "uniform"));
        const Vector3<real_t> initialVelocity = initParams.getParameter<Vector3<real_t>>("initialVelocity", Vector3<real_t>(0.0));

        const real_t kappa = parameters.getParameter<real_t>("kappa", real_t(0.42));
        const real_t B = parameters.getParameter<real_t>("B", real_t(5.5));
        const real_t roughnessLengthRatio = parameters.getParameter<real_t>("roughnessLengthRatio", real_t(1e-4));
        const real_t referenceHeight = parameters.getParameter<real_t>("referenceHeight_LU", real_t(-1));
        const uint32_t samplingHeight = parameters.getParameter<uint32_t>("samplingHeight_LU", uint32_t(0));
        const real_t roughnessLength = roughnessLengthRatio * referenceHeight;

        const uint_t welfordIntervalOutput = walberlaEnv.config()->getOneBlock("Output").getParameter<uint_t>("welfordIntervalOutput", uint_t(0));

        // u_tau from MOST
        const real_t u_tau = kappa * initialVelocity[0] / std::log(real_t(1) / roughnessLengthRatio);

        /// DOMAIN INITIALISATION

        WALBERLA_LOG_INFO_ON_ROOT("Initialise domain...")
        if(!domainSetup.snapshotHandler_->loadSnapshot) {

            std::unique_ptr<domain::DomainInitialisation> initialiser{};
            if (initType == domain::DomainInitialisation::UNIFORM) {
                initialiser = std::make_unique<domain::UniformInitialisation>(blocks, initialVelocity);
            } else if (initType == domain::DomainInitialisation::ASMUTH) {
                initialiser = std::make_unique<domain::AsmuthInitialisation>(blocks, u_tau, roughnessLength, kappa, B,
                                                                             viscosity_LU, domainSetup.domainSize_, true);
            } else if (initType == domain::DomainInitialisation::LOG_LAW) {
                // asmuth without perturbation velocities
                initialiser = std::make_unique<domain::AsmuthInitialisation>(blocks, u_tau, roughnessLength, kappa, B,
                                                                             viscosity_LU, domainSetup.domainSize_, false);
            }

            using MacroSetter_T = walberla::pystencils::waLBerlaWind_MacroSetter;
            MacroSetter_T setter(densityFieldID, forceFieldID, pdfFieldID, velocityFieldID);
            initialiser->setViaVelocityField<VectorField_T, MacroSetter_T>(velocityFieldID, setter);
            initialiser->setViaVelocityField<VectorField_T, MacroSetter_T>(meanVelocityFieldOutputID, setter);
            initialiser->setViaVelocityField<VectorField_T, MacroSetter_T>(meanVelocityFieldWFBID, setter);
//        initialiser->setPDFField<PdfField_T>(pdfFieldID);
        } else {
            if(domainSetup.snapshotHandler_->meanVelFile.empty()) {
                WALBERLA_ABORT("Mean velocity not loaded from snapshot but also not initialised.")
            }
            // initialise velocity field from pfd field if snapshot
            walberla::pystencils::waLBerlaWind_MacroGetter getter(densityFieldID, forceFieldID, pdfFieldID, velocityFieldID);
            for( auto & block : *blocks ) {
                getter(&block);
            }
        }

        SweepCollection_T sweepCollection( blocks, densityFieldID, forceFieldID, nutFieldID, pdfFieldID, velocityFieldID, omega );

        WALBERLA_MPI_BARRIER()
        WALBERLA_LOG_INFO_ON_ROOT("Initialisation done")

        /// BOUNDARIES

        WALBERLA_LOG_INFO_ON_ROOT("Setting up boundaries...")

        boundarySetup.fillFlagFieldFromConfig<FlagField_T>(blocks, flagFieldID, FluidFlagUID,
                                                           NoSlipFlagUID, WFBFlagUID, SymmetryFlagUID, UniformInflowFlagUID,
                                                           LogLawInflowFlagUID, OutflowFlagUID);

        if(boundarySetup.wallType() == WallSetup::WFB){
            WALBERLA_CHECK(referenceHeight > real_t(0), "Reference height must be given in the parameter file for wall function boundary conditions.")
            WALBERLA_CHECK(boundarySetup.environmentSetup() == EnvironmentSetup::Open,
                           "Wall-function bounce is currently only supported for an open environment")
            WALBERLA_CHECK(samplingHeight > uint_t(0), "Sampling height must be given in the parameter file for wall function boundary conditions.")
        }

        const auto& inflowVelocity = boundarySetup.inflowVelocity();

        std::function< walberla::math::Vector3< real_t >(const walberla::Cell&, const std::shared_ptr<walberla::StructuredBlockForest >&, walberla::IBlock&)>
            velocityInitialisation = std::bind(turbine_core::boundary::VelocityCallback, std::placeholders::_1,
                std::placeholders::_2, std::placeholders::_3, roughnessLength, kappa, u_tau);

        BoundaryCollection_T boundaryCollection( blocks, flagFieldID, pdfFieldID, FluidFlagUID, forceFieldID,
                                                 meanVelocityFieldWFBID, samplingHeight,
                                                 roughnessLengthRatio * referenceHeight,
                                                 inflowVelocity[0], inflowVelocity[1], inflowVelocity[2],
                                                 velocityInitialisation);

        const int periodicShiftValue = boundariesConfig.getParameter<int>("periodicShiftValue", int(real_t(0.33) * real_t(domainSetup.domainSize_[2])));
        walberla::boundary::ShiftedPeriodicity<PdfField_T> shiftedPeriodicity(blocks, pdfFieldID, fieldGhostLayers,
                                                                              0, 1, periodicShiftValue);

        /// set IDs for farm
        farm.setBlockForest(blocks);
        farm.setFieldIDs(densityFieldID, velocityFieldID, forceFieldID);


        /// WELFORD

        walberla::pystencils::waLBerlaWind_Welford welfordSweep(meanVelocityFieldWFBID, velocityFieldID, 0.0);
        walberla::pystencils::waLBerlaWind_WelfordSOP welfordSweepOutput(meanVelocityFieldOutputID, sumOfSquaresFieldID, velocityFieldID, 0.0);
        walberla::pystencils::waLBerlaWind_SoSResetter welfordSoSResetter(sumOfSquaresFieldID);

        /// TIMELOOP

        WALBERLA_LOG_INFO_ON_ROOT("Creating time loop...")

        // create time loop
        walberla::SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

        farm.setTimeloop(&timeloop);

        // driving forces
        walberla::wind::FlowDriverCollection flowDriver(blocks, &timeloop, globalConfig, domainSetup, forceFieldID, velocityFieldID, fieldGhostLayers);
        for( auto & block : *blocks) {
            flowDriver(&block);
        }

        // create communication for fields
        walberla::blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication( blocks );
        communication.addPackInfo( std::make_shared< walberla::lbm_generated::UniformGeneratedPdfPackInfo<PdfField_T> >( pdfFieldID ) );
        communication.addPackInfo( std::make_shared< walberla::field::communication::PackInfo<ScalarField_T> >( densityFieldID ) );
        communication.addPackInfo( std::make_shared< walberla::field::communication::PackInfo<VectorField_T > >( velocityFieldID ) );

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
        auto applyControl = [&farm](){farm.applyControl();};
        auto syncForces = [&](){farm.syncNextNeighbour(mpi::SYNC_FORCE, maxKernelWidth+uint_t(real_t(diameter)/real_t(2.0)), false);};
        auto syncMacros = [&](){farm.syncNextNeighbour(mpi::SYNC_MACRO, 1, false);};
        auto spreadForces = [&farm](walberla::IBlock * block){farm.spreadForces(block);};

        // initial communication
        communication();

        // add LBM sweep and communication to time loop
        // output
        TurbineOutput_T turbineOutput{
                blocks, &timeloop,
                globalConfig->getOneBlock("Output"), FluidFlagUID,
                fieldMap
        };

        //FIXME the current implementation expects the omega field, not the nu_t field. Strain rate writer must be called BEFORE eddy viscosity writer

        walberla::pystencils::waLBerlaWind_StrainRateWriter strainRateWriter(forceFieldID, nutFieldID, pdfFieldID, strainRateFieldID);
        turbineOutput.addBeforeFunction( [&]() {
            for( auto & block : *blocks ) {
                strainRateWriter(&block);
            }});

        walberla::pystencils::waLBerlaWind_EddyViscosityWriter eddyViscosityWriter(nutFieldID, omega);
        turbineOutput.addBeforeFunction( [&]() {
            for( auto & block : *blocks ) {
                eddyViscosityWriter(&block);
            }});

        auto output = std::bind( &TurbineOutput_T::write, &turbineOutput );

        //NOTE must convert sweeps that are altered to lambdas, otherwise copy and counter will stay 0
        welfordSweep.setCounter(real_t(initialWelfordCounter));
        welfordSweepOutput.setCounter(real_t(initialWelfordCounter));
        auto welfordLambda = [&welfordSweep](walberla::IBlock * block) {
            welfordSweep(block);
        };
        auto welfordLambdaOutput = [&welfordSweepOutput](walberla::IBlock * block) {
            welfordSweepOutput(block);
        };

        timeloop.add() << walberla::BeforeFunction( communication, "field communication" )
                       << walberla::BeforeFunction( [&boundarySetup, &shiftedPeriodicity]() {
                            if( boundarySetup.inflowType() == InflowSetup::ShiftedPeriodic )
                                shiftedPeriodicity();
                       }, "Shifted periodicity" )
                       << walberla::Sweep( boundaryCollection.getSweep(), "Boundary handling" );
        timeloop.add() << walberla::Sweep( sweepCollection.stream() , "LBM stream sweep");
        timeloop.add() << walberla::Sweep( sweepCollection.collide() , "LBM collide sweep")
                       << walberla::AfterFunction( output, "VTK output" )
                       << walberla::AfterFunction( [&](){
                           domainSetup.snapshotHandler_->createSnapshot(timeloop.getCurrentTimeStep(), blocks, pdfFieldID, forceFieldID,
                                                                        uint_t(welfordSweep.getCounter()), meanVelocityFieldWFBID,
                                                                        sumOfSquaresFieldID);
                       }, "Snapshot creation");
        if(boundarySetup.wallType() == WallSetup::WFB) {
            timeloop.add() << walberla::BeforeFunction([&](){
                                  welfordSweep.setCounter(real_t(welfordSweep.getCounter()+1));
                              }, "welford sweep")
                           << walberla::Sweep(welfordLambda, "welford sweep");
        }

        timeloop.add() << walberla::BeforeFunction([&](){
                              if(welfordIntervalOutput // && int(welfordSweepOutput.counter_)
                                 && (uint_t(welfordSweepOutput.getCounter()) % welfordIntervalOutput == 0)){
                                  welfordSweepOutput.setCounter(real_t(1));
                                  // reset the field: swap pointers

                                  // Take care this is called only after the VTK output, otherwise we write the wrong data in the meanVelocity output
                                  // We live a dangerous life: we swap the velocity field counter
                                  for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) {
                                      auto velocityField = blockIt->getData<VectorField_T>( velocityFieldID );
                                      auto meanVelocityField = blockIt->getData<VectorField_T>( meanVelocityFieldOutputID );
                                      velocityField->swapDataPointers(meanVelocityField);

                                      welfordSoSResetter(blockIt.get());
                                  }
                              } else {
                                  welfordSweepOutput.setCounter(real_t(welfordSweepOutput.getCounter()+1));
                              }
                          }, "welford sweep Output")
                       << walberla::Sweep(welfordLambdaOutput, "welford sweep Output");

        timeloop.add() << walberla::BeforeFunction( [&]() {
                              farm.callback(Component::Function::UPDATE_DISCRETISATION);
                          }, "Rotation" )
                       << walberla::Sweep( evaluateDensityAndVelocity, "Evaluate density and velocity", sets::turbineSelector() )
                       << walberla::AfterFunction(syncMacros, "Communicate macroscopic turbine data")
                       << walberla::AfterFunction( applyControl , "Apply Control");

        timeloop.add() << walberla::Sweep( flowDriver, "Setting driving force");
        timeloop.add() << walberla::BeforeFunction( [&]() {
                              farm.callback( Component::Function::CALCULATE_FORCES );
                          }, "Calculate forces" )
                       << walberla::BeforeFunction( syncForces, "Communicate force data" )
                       << walberla::Sweep(spreadForces, "Spread forces", sets::turbineSelector() )
                       << walberla::AfterFunction( [&]() {
                           farm.callback(Component::Function::RESET_COMMUNICATION_DATA);
                       }, "ResetCommunicationData" );

        // LBM stability check
        timeloop.addFuncAfterTimeStep( walberla::makeSharedFunctor( walberla::field::makeStabilityChecker< PdfField_T, FlagField_T >(globalConfig, blocks, pdfFieldID, flagFieldID, FluidFlagUID ) ),
                                       "LBM stability check" );

        // log remaining time
        timeloop.addFuncAfterTimeStep( walberla::timing::RemainingTimeLogger( timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency ), "remaining time logger" );

        WALBERLA_LOG_INFO_ON_ROOT("Running timeloop...")

        WALBERLA_MPI_WORLD_BARRIER()

        walberla::WcTimingPool timing;
        timing.registerTimer("Force Output");

        walberla::WcTimer timer;
        timer.start();

        for( uint_t i = 0; i < timesteps; ++i ) {

            timeloop.singleStep(timing);

            timing["Force Output"].start();
            // farm.callback(Component::Output::ORIENTATIONS);
            for(auto & block : *blocks) {
                farm.writeForceOutput(&block, forceOutputStart, forceOutputFrequency);
            }
            farm.writeTurbinePowerAndThrust();
            timing["Force Output"].end();
        }

        timer.end();

        double time = timer.max();
        walberla::mpi::reduceInplace(time, walberla::mpi::MAX);

        const auto timeloopTiming = timing.getReduced();
        WALBERLA_LOG_INFO_ON_ROOT("Timeloop timing:\n" << *timeloopTiming)

        walberla::lbm::PerformanceEvaluation<FlagField_T> performance(blocks, flagFieldID, FluidFlagUID);
        performance.logResultOnRoot(timesteps, time);

        farm.calculatePowerAndThrust();

        logging::WindTurbineApplicationLogger<domain::DomainSetup, domain::BoundarySetup, codegen::KernelInfo, walberla::timing::WcPolicy> logger {
            domainSetup, boundarySetup, windFarmConfig, &timing, performance.loggingString(timesteps, time)
        };

        return EXIT_SUCCESS;

    }

}

int main(int argc, char** argv) {
    return turbine_core::main(argc, argv);
}
