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

    static const uint_t fieldGhostLayers =  2 ;

    // define different interpolators and distributors

    template<class Kernel_T>
    using KernelFieldInterpolator_T = projectors::KernelFieldInterpolator<field::Field<real_t>, Kernel_T>;

    template<class Kernel_T>
    using KernelDistributor_T = projectors::KernelFieldDistributor<field::Field<real_t>, Kernel_T>;

    using InterpolationKernel_T = projectors::SmoothedDiracDeltaKernel;
    using DistributionKernel_T = projectors::SmoothedDiracDeltaKernel;

    using WindFarm_T = WindFarm<topology::cpu::Tree, creator::ConfigTurbineCreator,
            ScalarField_T, VectorField_T, VectorField_T,
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

//        projectors::GaussianFunction5<>::setSigma(real_t(7.07106781));

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
        
        WALBERLA_ASSERT(!domainSetup.snapshotHandler_->loadSnapshot, "Restart is currently not supported for non-uniform grids.")

        refinement::StaticRefinementHandler<WindFarm_T> staticRefinementHandler{ globalConfig->getOneBlock("DomainSetup"), &farm, boundarySetup.environmentSetup() };
        auto blocks = domainSetup.createStructuredBlockForest(staticRefinementHandler, fieldGhostLayers);

        farm.setBlockForest(blocks);
        farm.initialiseForceAndControlModel();

        /// FIELDS

        WALBERLA_LOG_INFO_ON_ROOT("Creating fields...")

        //layout - DO NOT CHANGE
        auto layout = codegen::KernelInfo::layout;
        const StorageSpecification_T storageSpec{};

        std::map<output::Fields::Types, BlockDataID> fieldMap {};
        uint_t initialWelfordCounter = 0;

        auto allocator = std::make_shared< walberla::field::AllocateAligned<real_t, 64> >();
        BlockDataID forceFieldCpuID;
        BlockDataID meanVelocityOutputFieldCpuID;
        BlockDataID meanVelocityWfbFieldCpuID;
        BlockDataID pdfFieldCpuID;
        BlockDataID sumOfCubesFieldCpuID;
        BlockDataID sumOfSquaresFieldCpuID;

        if(domainSetup.snapshotHandler_->loadSnapshot) {

        	auto & ssh = domainSetup.snapshotHandler_;
        	
        	auto vectorFieldDataHandling = std::make_shared<field::FlattenedFieldHandling<VectorField_T>> (blocks, fieldGhostLayers, layout, allocator);
        	const auto forceFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->forceFile).c_str());
        	
        	forceFieldCpuID = blocks->loadBlockData(forceFileName, vectorFieldDataHandling, "force field CPU");
        	
        	auto pdfFieldDataHandling = std::make_shared< walberla::lbm_generated::internal::PdfFieldHandling< StorageSpecification_T > >( blocks, storageSpec, fieldGhostLayers, layout, allocator );
        	
        	const auto pdfFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->pdfFile).c_str());
        	pdfFieldCpuID = blocks->loadBlockData(pdfFileName, pdfFieldDataHandling, "pdf field CPU");
        	

        	if(!ssh->meanVelFile.empty()) {
        		const auto meanVelFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->meanVelFile).c_str());
        		initialWelfordCounter = ssh->welfordCounter;

        		meanVelocityWfbFieldCpuID    = blocks->loadBlockData(meanVelFileName, vectorFieldDataHandling, "mean velocity wfb field CPU");
        		meanVelocityOutputFieldCpuID = blocks->loadBlockData(meanVelFileName, vectorFieldDataHandling, "mean velocity output field CPU");

        		if(!ssh->sosFile.empty()) {
        			auto secondOrderTensorFieldDataHandling = std::make_shared<field::FlattenedFieldHandling<SecondOrderTensorField_T>> (blocks, fieldGhostLayers, layout, allocator);
        			const auto sosFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->sosFile).c_str());
        			sumOfSquaresFieldCpuID = blocks->loadBlockData(sosFileName, secondOrderTensorFieldDataHandling, "sum of squares field");
        		}
        		if(!ssh->socFile.empty()) {
        			auto thirdOrderTensorFieldDataHandling = std::make_shared<field::FlattenedFieldHandling<ThirdOrderTensorField_T>> (blocks, fieldGhostLayers, layout, allocator);
        			const auto socFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->socFile).c_str());
        			sumOfCubesFieldCpuID = blocks->loadBlockData(socFileName, thirdOrderTensorFieldDataHandling, "sum of cubes field CPU");
        		}
        	}
        } else {
        	forceFieldCpuID = walberla::field::addToStorage<VectorField_T>(blocks, "force field CPU", real_t(0), layout, fieldGhostLayers, allocator);
        	meanVelocityOutputFieldCpuID = walberla::field::addToStorage<VectorField_T>(blocks, "mean velocity output field CPU", real_t(0), layout, fieldGhostLayers, allocator);
        	meanVelocityWfbFieldCpuID = walberla::field::addToStorage<VectorField_T>(blocks, "mean velocity wfb field CPU", real_t(0), layout, fieldGhostLayers, allocator);
        	pdfFieldCpuID = walberla::lbm_generated::addPdfFieldToStorage(blocks, "pdf field CPU", storageSpec, fieldGhostLayers, layout, walberla::Set<walberla::SUID>::emptySet(), walberla::Set<walberla::SUID>::emptySet(), allocator);
        	sumOfCubesFieldCpuID = walberla::field::addToStorage<ThirdOrderTensorField_T>(blocks, "sum of cubes field CPU", real_t(0), layout, fieldGhostLayers, allocator);
        	sumOfSquaresFieldCpuID = walberla::field::addToStorage<SecondOrderTensorField_T>(blocks, "sum of squares field CPU", real_t(0), layout, fieldGhostLayers, allocator);
        }

        BlockDataID densityFieldCpuID = walberla::field::addToStorage<ScalarField_T>(blocks, "density field CPU", real_t(1), layout, fieldGhostLayers, allocator);
        BlockDataID eddyViscosityFieldCpuID = walberla::field::addToStorage<ScalarField_T>(blocks, "eddy viscosity field CPU", real_t(0), layout, fieldGhostLayers, allocator);
        BlockDataID flagFieldID = walberla::field::addFlagFieldToStorage<FlagField_T>(blocks, "flag field", fieldGhostLayers);
        BlockDataID meanEddyViscosityFieldCpuID = walberla::field::addToStorage<ScalarField_T>(blocks, "mean eddy viscosity field CPU", real_t(0), layout, fieldGhostLayers, allocator);
        BlockDataID meanStrainRateFieldCpuID = walberla::field::addToStorage<SecondOrderTensorField_T>(blocks, "mean strain rate field CPU", real_t(0), layout, fieldGhostLayers, allocator);
        BlockDataID omegaFieldCpuID = walberla::field::addToStorage<ScalarField_T>(blocks, "omega field CPU", real_t(0), layout, fieldGhostLayers, allocator);
        BlockDataID strainRateFieldCpuID = walberla::field::addToStorage<SecondOrderTensorField_T>(blocks, "strain rate field CPU", real_t(0), layout, fieldGhostLayers, allocator);
        BlockDataID velocityFieldCpuID = walberla::field::addToStorage<VectorField_T>(blocks, "velocity field CPU", real_t(0), layout, fieldGhostLayers, allocator);


        fieldMap[output::Fields::DENSITY] = densityFieldCpuID;
        fieldMap[output::Fields::EDDY_VISCOSITY] = eddyViscosityFieldCpuID;
        fieldMap[output::Fields::FLAG] = flagFieldID;
        fieldMap[output::Fields::FORCE] = forceFieldCpuID;
        fieldMap[output::Fields::MEAN_EDDY_VISCOSITY] = meanEddyViscosityFieldCpuID;
        fieldMap[output::Fields::MEAN_STRAIN_RATE] = meanStrainRateFieldCpuID;
        fieldMap[output::Fields::MEAN_VELOCITY_OUTPUT] = meanVelocityOutputFieldCpuID;
        fieldMap[output::Fields::MEAN_VELOCITY_WFB] = meanVelocityWfbFieldCpuID;
        fieldMap[output::Fields::OMEGA] = omegaFieldCpuID;
        fieldMap[output::Fields::PDF] = pdfFieldCpuID;
        fieldMap[output::Fields::STRAIN_RATE] = strainRateFieldCpuID;
        fieldMap[output::Fields::SUM_OF_CUBES] = sumOfCubesFieldCpuID;
        fieldMap[output::Fields::SUM_OF_SQUARES] = sumOfSquaresFieldCpuID;
        fieldMap[output::Fields::VELOCITY] = velocityFieldCpuID;

        /// ABL calculations

        auto initParams = globalConfig->getOneBlock("Initialisation");
        const auto initType = domain::DomainInitialisation::toType(initParams.getParameter<std::string>("type"));
        const Vector3<real_t> initialVelocity = initParams.getParameter<Vector3<real_t>>("initialVelocity");

        const real_t kappa = parameters.getParameter<real_t>("kappa", real_t(0.42));
        const real_t B = parameters.getParameter<real_t>("B", real_t(5.5));
        const real_t roughnessLengthRatio = parameters.getParameter<real_t>("roughnessLengthRatio", real_t(1e-4));
        const real_t referenceHeight = parameters.getParameter<real_t>("referenceHeight_LU", real_t(-1));
        const uint32_t samplingHeight = parameters.getParameter<uint32_t>("samplingHeight_LU", uint32_t(0));
        const real_t roughnessLength = roughnessLengthRatio * referenceHeight;

        const uint_t welfordInterval = walberlaEnv.config()->getOneBlock("Output").getParameter<uint_t>("welfordInterval", uint_t(0));
        // uTau from MOST
        const real_t uTau = kappa * initialVelocity[0] / std::log(real_t(1) / roughnessLengthRatio);

        /// DOMAIN INITIALISATION

        WALBERLA_LOG_INFO_ON_ROOT("Initialise domain...")
        if(!domainSetup.snapshotHandler_->loadSnapshot) {
        			
            std::unique_ptr<domain::DomainInitialisation> initialiser{};
            if (initType == domain::DomainInitialisation::UNIFORM) {
                initialiser = std::make_unique<domain::UniformInitialisation>(blocks, initialVelocity);
            } else if (initType == domain::DomainInitialisation::ASMUTH) {
                initialiser = std::make_unique<domain::AsmuthInitialisation>(blocks, uTau, roughnessLength, kappa, B,
                                                                             viscosity_LU, domainSetup.domainSize_, true);
            } else if (initType == domain::DomainInitialisation::LOG_LAW) {
                // asmuth without perturbation velocities
                initialiser = std::make_unique<domain::AsmuthInitialisation>(blocks, uTau, roughnessLength, kappa, B,
                                                                             viscosity_LU, domainSetup.domainSize_, false);
            }
            using PdfSetter_T = walberla::pystencils::waLBerlaWind_PdfSetter;
            PdfSetter_T setter(densityFieldCpuID, forceFieldCpuID, pdfFieldCpuID, velocityFieldCpuID);
            initialiser->setViaVelocityField<VectorField_T, PdfSetter_T>(velocityFieldCpuID, setter);

        	initialiser->setViaVelocityField<VectorField_T, PdfSetter_T>(meanVelocityOutputFieldCpuID, setter);
        	initialiser->setViaVelocityField<VectorField_T, PdfSetter_T>(meanVelocityWfbFieldCpuID, setter);
        } else {

        	if(domainSetup.snapshotHandler_->meanVelFile.empty()) {
        		WALBERLA_ABORT("Mean velocity not loaded from snapshot but also not initialised.")
        	}

        	// initialise velocity field from pfd field if snapshot
        	walberla::pystencils::waLBerlaWind_MacroGetter getter(densityFieldCpuID, forceFieldCpuID, pdfFieldCpuID, velocityFieldCpuID);
        	    for( auto & block : *blocks ) {
                            getter(&block);
                        }
        }

        SweepCollection_T sweepCollection( blocks, densityFieldCpuID, eddyViscosityFieldCpuID, forceFieldCpuID, omegaFieldCpuID, 
                                           pdfFieldCpuID, velocityFieldCpuID, omega );

        // for (auto& block : *blocks) {
        //     sweepCollection.initialise(&block);
        //     sweepCollection.initialise(&block, fieldGhostLayers);
        // }
        WALBERLA_MPI_BARRIER()
        WALBERLA_LOG_INFO_ON_ROOT("Initialisation done")

        /// BOUNDARIES

        WALBERLA_LOG_INFO_ON_ROOT("Setting up boundaries...")

        boundarySetup.fillFlagFieldFromConfig<FlagField_T>(blocks, flagFieldID, FluidFlagUID,
                                                           NoSlipFlagUID, WFBFlagUID, SymmetryFlagUID, UniformInflowFlagUID,
                                                           LogLawInflowFlagUID, OutflowFlagUID);

        if(boundarySetup.wallType() == WallSetup::WFB){
            WALBERLA_CHECK(referenceHeight > real_t(0), "Reference height must be given in the parameter file for wall function boundary conditions.")
            WALBERLA_CHECK(boundarySetup.environmentSetup() != EnvironmentSetup::Tunnel,
                           "Wall-function bounce is currently not supported for the tunnel environment")
            WALBERLA_CHECK(samplingHeight > uint_t(0), "Sampling height must be given in the parameter file for wall function boundary conditions.")
        }
        

        const auto& inflowVelocity = boundarySetup.inflowVelocity();
        auto velocityInit = boundary::velocityInit(roughnessLength, kappa, uTau);
        BoundaryCollection_T boundaryCollection( 
        	blocks, flagFieldID, pdfFieldCpuID, FluidFlagUID, forceFieldCpuID, 
        	meanVelocityWfbFieldCpuID, samplingHeight,
        	roughnessLengthRatio * referenceHeight,
        	inflowVelocity[0], inflowVelocity[1], inflowVelocity[2], velocityInit
        );

        const int periodicShiftValue = boundariesConfig.getParameter<int>("periodicShiftValue", int(real_t(0.33) * real_t(domainSetup.domainSize_[2])));
        walberla::boundary::ShiftedPeriodicity<PdfField_T> shiftedPeriodicity (blocks, pdfFieldCpuID, fieldGhostLayers, 0, 1, periodicShiftValue);

        /// set IDs for farm
        farm.setFieldIDs(densityFieldCpuID, velocityFieldCpuID, forceFieldCpuID);

        /// field resetter
        
        walberla::pystencils::waLBerlaWind_SoSResetter welfordOutputSosResetter(sumOfSquaresFieldCpuID);

        walberla::pystencils::waLBerlaWind_SoCResetter welfordOutputSocResetter(sumOfCubesFieldCpuID);

        /// WELFORDWFB
        walberla::pystencils::waLBerlaWind_WelfordWFB welfordWFBSweep(
        	meanVelocityWfbFieldCpuID, velocityFieldCpuID, 
        	real_t(initialWelfordCounter));
        /// WELFORDOUTPUT
        walberla::pystencils::waLBerlaWind_WelfordOutput welfordOutputSweep(
        	meanVelocityOutputFieldCpuID, sumOfCubesFieldCpuID, sumOfSquaresFieldCpuID, velocityFieldCpuID, 
        	real_t(initialWelfordCounter));
        /// WELFORDEDDYVISCOSITY
        walberla::pystencils::waLBerlaWind_WelfordEddyViscosity welfordEddyViscositySweep(
        	eddyViscosityFieldCpuID, meanEddyViscosityFieldCpuID, 
        	real_t(initialWelfordCounter));
        /// WELFORDSTRAINRATE
        walberla::pystencils::waLBerlaWind_WelfordStrainRate welfordStrainRateSweep(
        	meanStrainRateFieldCpuID, strainRateFieldCpuID, 
        	real_t(initialWelfordCounter));

        /// COMMUNICATION

        WALBERLA_LOG_INFO_ON_ROOT("Set up communication...")

        // create communication for fields
        auto communication = std::make_shared<walberla::blockforest::communication::NonUniformBufferedScheme<CommunicationStencil_T>>(blocks);
        auto pdfPackInfo = walberla::lbm_generated::setupNonuniformPdfCommunication<PdfField_T>(blocks, pdfFieldCpuID);
        communication->addPackInfo( pdfPackInfo );
        auto densityPackInfo = std::make_shared< walberla::field::refinement::PackInfo<ScalarField_T, CommunicationStencil_T> >(densityFieldCpuID);
        communication->addPackInfo( densityPackInfo );
        auto velocityPackInfo = std::make_shared< walberla::field::refinement::PackInfo<VectorField_T, CommunicationStencil_T> >(velocityFieldCpuID);
        communication->addPackInfo( velocityPackInfo );

        /// TIMELOOP

        WALBERLA_LOG_INFO_ON_ROOT("Creating time loop...")

        using TimingPolicy_T = walberla::timing::WcPolicy;
        // create time loop
        auto timeloop = walberla::timeloop::SweepTimeloop<TimingPolicy_T>(blocks->getBlockStorage(), timesteps);

        farm.setTimeloop(&timeloop);

        /// driving forces
        walberla::wind::FlowDriverCollection flowDriver(blocks, &timeloop, globalConfig, domainSetup, forceFieldCpuID, velocityFieldCpuID, fieldGhostLayers);
        for(auto & block : *blocks) { flowDriver(&block); }

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

        // add LBM sweep and communication to time loop
        // output
        TurbineOutput_T turbineOutput{
                blocks, &timeloop,
                globalConfig->getOneBlock("Output"), FluidFlagUID,
                fieldMap
        };

        walberla::pystencils::waLBerlaWind_StrainRateWriter strainRateWriter(forceFieldCpuID, omegaFieldCpuID, pdfFieldCpuID, strainRateFieldCpuID);

        turbineOutput.addBeforeFunction( [&]() {
        for( auto & block : *blocks ) {
            strainRateWriter(&block);
        }});

        auto output = [&turbineOutput]() { turbineOutput.write(); };
        

        auto welfordWFBLambda = [&welfordWFBSweep](walberla::IBlock * block) {
            welfordWFBSweep(block);
        };
        auto welfordOutputLambda = [&welfordOutputSweep](walberla::IBlock * block) {
            welfordOutputSweep(block);
        };
        auto welfordEddyViscosityLambda = [&welfordEddyViscositySweep](walberla::IBlock * block) {
            welfordEddyViscositySweep(block);
        };
        auto welfordStrainRateLambda = [&welfordStrainRateSweep](walberla::IBlock * block) {
            welfordStrainRateSweep(block);
        };

        refinement::CustomRecursiveTimeStep<PdfField_T, SweepCollection_T, BoundaryCollection_T> recursiveTimestep (
        	blocks, pdfFieldCpuID, sweepCollection, boundaryCollection, communication, pdfPackInfo );

        recursiveTimestep.addAfterCollideFunction(output, "VTK output", 0);
        recursiveTimestep.addAfterCollideFunction( [&](){
                          if ((domainSetup.snapshotHandler_->storeSnapshot) && (timeloop.getCurrentTimeStep() % domainSetup.snapshotHandler_->snapshotFrequency == 0) && (timeloop.getCurrentTimeStep() > 0)) {
                             domainSetup.snapshotHandler_->createSnapshot(timeloop.getCurrentTimeStep(), blocks, pdfFieldCpuID, forceFieldCpuID,
                                                                          uint_t(welfordWFBSweep.getCounter()), meanVelocityWfbFieldCpuID, sumOfSquaresFieldCpuID, sumOfCubesFieldCpuID);
                          }
                          }, "Snapshot creation");
        if(boundarySetup.wallType() == WallSetup::WFB) {
        	recursiveTimestep.addAfterCollideFunction([&]() { welfordWFBSweep.setCounter(real_t(welfordWFBSweep.getCounter() + 1)); }, "Welford Sweep");
        	recursiveTimestep.addAfterCollideFunction([&]() {
        for (auto &block: *blocks) {
        	welfordWFBLambda(&block);
        }}, "Welford Sweep");

        }
        	recursiveTimestep.addAfterCollideFunction([&](){
                           	if(welfordInterval && (uint_t(welfordOutputSweep.getCounter()) % welfordInterval == 0)){
                           		welfordOutputSweep.setCounter(real_t(1));
                           		for ( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) {
                           			auto dst = blockIt->getData<VectorField_T>( meanVelocityOutputFieldCpuID );
                           			const auto src = blockIt->getData<VectorField_T>( velocityFieldCpuID );
                           			field::fieldCopy( dst, src );
                           			welfordOutputSosResetter(blockIt.get());
                           			welfordOutputSocResetter(blockIt.get());
                           		}
                           	} else {
                           		welfordOutputSweep.setCounter(real_t(welfordOutputSweep.getCounter() + 1));
                           	}
                           }, "Welford Output Sweep");
        	recursiveTimestep.addAfterCollideFunction([&]() {
        for (auto &block: *blocks) {
        	welfordOutputLambda(&block);
        }}, "Welford Output Sweep");

        	recursiveTimestep.addAfterCollideFunction([&](){
                           	if(welfordInterval && (uint_t(welfordEddyViscositySweep.getCounter()) % welfordInterval == 0)){
                           		welfordEddyViscositySweep.setCounter(real_t(1));
                           		for ( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) {
                           			auto dst = blockIt->getData<ScalarField_T>( meanEddyViscosityFieldCpuID );
                           			const auto src = blockIt->getData<ScalarField_T>( eddyViscosityFieldCpuID );
                           			field::fieldCopy( dst, src );
                           		}
                           	} else {
                           		welfordEddyViscositySweep.setCounter(real_t(welfordEddyViscositySweep.getCounter() + 1));
                           	}
                           }, "Welford Eddy Viscosity Sweep");
        	recursiveTimestep.addAfterCollideFunction([&]() {
        for (auto &block: *blocks) {
        	welfordEddyViscosityLambda(&block);
        }}, "Welford Eddy Viscosity Sweep");

        	recursiveTimestep.addAfterCollideFunction([&](){
                           	if(welfordInterval && (uint_t(welfordStrainRateSweep.getCounter()) % welfordInterval == 0)){
                           		welfordStrainRateSweep.setCounter(real_t(1));
                           		for ( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) {
                           			auto dst = blockIt->getData<SecondOrderTensorField_T>( meanStrainRateFieldCpuID );
                           			const auto src = blockIt->getData<SecondOrderTensorField_T>( strainRateFieldCpuID );
                           			field::fieldCopy( dst, src );
                           		}
                           	} else {
                           		welfordStrainRateSweep.setCounter(real_t(welfordStrainRateSweep.getCounter() + 1));
                           	}
                           }, "Welford Strain Rate Sweep");
        	recursiveTimestep.addAfterCollideFunction([&]() {
        for (auto &block: *blocks) {
        	welfordStrainRateLambda(&block);
        }}, "Welford Strain Rate Sweep");

        recursiveTimestep.addAfterCollideFunction([&]() { farm.callback(Component::Function::UPDATE_DISCRETISATION); }, "Rotation");
        recursiveTimestep.addAfterCollideFunction([&]() {
        for (auto &block: *blocks) {
        	evaluateDensityAndVelocity(&block);
        }}, "Evaluate density and velocity");

        recursiveTimestep.addAfterCollideFunction(syncMacros, "Communicate macroscopic turbine data");
        recursiveTimestep.addAfterCollideFunction(applyControl, "Apply Control");
        recursiveTimestep.addAfterCollideFunction([&]() {
        for (auto &block: *blocks) {
        	flowDriver(&block);
        }}, "Setting driving force");

        recursiveTimestep.addAfterCollideFunction([&]() { farm.callback( Component::Function::CALCULATE_FORCES ); }, "Calculate forces");
        recursiveTimestep.addAfterCollideFunction(syncForces, "Communicate force data");
        recursiveTimestep.addAfterCollideFunction([&]() {
        for (auto &block: *blocks) {
        	spreadForces(&block);
        }}, "Spread forces");

        recursiveTimestep.addAfterCollideFunction([&]() { farm.callback(Component::Function::RESET_COMMUNICATION_DATA); }, "Reset communication data");

        recursiveTimestep.addRefinementToTimeLoop( timeloop );

        // LBM stability check
        timeloop.addFuncAfterTimeStep( walberla::makeSharedFunctor( walberla::field::makeStabilityChecker< PdfField_T, FlagField_T >(globalConfig, blocks, pdfFieldCpuID, flagFieldID, FluidFlagUID )), "LBM stability check" );

        // log remaining time
        timeloop.addFuncAfterTimeStep( walberla::timing::RemainingTimeLogger( timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency ), "Remaining time logger" );

        WALBERLA_LOG_INFO_ON_ROOT("Running timeloop...")

        WALBERLA_MPI_WORLD_BARRIER()

        walberla::timing::TimingPool<TimingPolicy_T> timing;
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

        logging::WindTurbineApplicationLogger<domain::DomainSetup, domain::BoundarySetup, codegen::KernelInfo, TimingPolicy_T> logger {
            domainSetup, boundarySetup, windFarmConfig, &timing, performance.loggingString(timesteps, time)
        };

        return EXIT_SUCCESS;

    } // main

} // namespace turbine_core

int main(int argc, char** argv) {
    return turbine_core::main(argc, argv);
}