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
#include "walberla_helper/field/FieldHandling.h"

#include "domain/BoundarySetup.h"
#include "domain/BoundarySetter.h"
#include "domain/BoundaryHandling.h"
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
#include "boundary/ShiftedPeriodicity.h"

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


    void copyBoundaryData( real_t * precursorField, real_t * concurrentField,
        std::vector<uint_t> cellIndices, const uint_t nbCells, const uint_t q, const uint_t fStride ) {
        for (uint_t id = 0; id < nbCells; ++id){
            for (uint_t qIdx = 0; qIdx < q; ++qIdx) {
                concurrentField[cellIndices[id] + qIdx*fStride] = precursorField[cellIndices[id] + qIdx*fStride];
            }
        }
    }

    int main(int argc, char** argv) {

        walberla::Environment walberlaEnv(argc, argv);

        walberla::logging::Logging::instance()->setLogLevel(walberla::logging::Logging::INFO);

        /// general simulation parameters
        auto globalConfig = walberlaEnv.config();
        auto parameters = globalConfig->getOneBlock("Parameters");

        uint_t timesteps = parameters.getParameter<uint_t>("timesteps", uint_t(10)); ++timesteps;
        const uint_t concurrentStart = parameters.getParameter< uint_t >("concurrentStart", uint_t(5000));
        const real_t remainingTimeLoggerFrequency = parameters.getParameter<real_t>("remainingTimeLoggerFrequency", real_t(3.0)); // in seconds

        const real_t length_SI = parameters.getParameter<real_t>("length_SI");
        const real_t length_LU = parameters.getParameter<real_t>("length_LU");
        const real_t velocity_SI = parameters.getParameter<real_t>("velocity_SI");
        const real_t velocity_LU = parameters.getParameter<real_t>("velocity_LU");
        const real_t viscosity_SI = parameters.getParameter<real_t>("viscosity_SI");
        const real_t density_SI = parameters.getParameter<real_t>("density_SI");

        Conversion::calculateConversionFactors(length_SI, length_LU, velocity_SI, velocity_LU, density_SI);

        Conversion::print();
        /// parameter conversion
        real_t viscosity_LU = (viscosity_SI / density_SI) * Conversion::C_t() / Conversion::C_l() / Conversion::C_l();

        WALBERLA_CHECK_GREATER(viscosity_LU, real_t(0.0), "Negative LU velocity - Did you calculate the conversion factors?")

        projectors::GaussianFunction5<>::setSigma(real_t(7.07106781));

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

        auto precursorBoundariesConfig = globalConfig->getOneBlock( "PrecursorBoundaries" );
        domain::BoundarySetup precursorBoundarySetup(precursorBoundariesConfig);
        domain::DomainSetup precursorDomainSetup ( globalConfig, precursorBoundarySetup.periodicity() );
        auto precursorBlocks = precursorDomainSetup.createUniformBlockForest(farm.getTurbineAABBs());


        auto concurrentBoundariesConfig = globalConfig->getOneBlock( "ConcurrentBoundaries" );
        domain::BoundarySetup concurrentBoundarySetup(concurrentBoundariesConfig);
        domain::DomainSetup concurrentDomainSetup ( globalConfig, concurrentBoundarySetup.periodicity() );
        auto concurrentBlocks = concurrentDomainSetup.createUniformBlockForest(farm.getTurbineAABBs());

        farm.setBlockForest(concurrentBlocks);
        farm.initialiseForceAndControlModel();

        /// FIELDS

        WALBERLA_LOG_INFO_ON_ROOT("Creating fields...")

        std::map<output::Fields::Types, BlockDataID> concurrentFieldMap {};
        std::map<output::Fields::Types, BlockDataID> precursorFieldMap {};

        //layout - DO NOT CHANGE
        auto layout = codegen::KernelInfo::layout;

        BlockDataID precursorPdfFieldID;
        BlockDataID precursorForceFieldID;
        BlockDataID precursorMeanVelocityFieldWFBID;
        BlockDataID precursorDensityFieldID = walberla::field::addToStorage<ScalarField_T>(precursorBlocks, "density field", real_t(1), layout, fieldGhostLayers);
        BlockDataID precursorVelocityFieldID = walberla::field::addToStorage<VectorField_T>(precursorBlocks, "precursor velocity field", real_t(0), layout, fieldGhostLayers);
        BlockDataID precursorFlagFieldID = walberla::field::addFlagFieldToStorage<FlagField_T>(precursorBlocks, "precursor flag field", fieldGhostLayers);
        BlockDataID precursorDummyNutFieldID = walberla::field::addToStorage<ScalarField_T>(precursorBlocks, "precursor eddy viscosity field", real_t(0), layout, fieldGhostLayers);

        // WARNING: the concurrent domain does not get an additional meanVelcityFieldWBF but uses the one from the precursor.
        // This needs to be kept in mind
        BlockDataID concurrentPdfFieldID;
        BlockDataID concurrentForceFieldID;
        BlockDataID concurrentMeanVelocityFieldOutputID;
        BlockDataID concurrentSumOfSquaresFieldID;
        BlockDataID concurrentMeanVelocityFieldWFBID;
        BlockDataID concurrentDensityFieldID = walberla::field::addToStorage<ScalarField_T>(concurrentBlocks, "density field", real_t(1), layout, fieldGhostLayers);
        BlockDataID concurrentNutFieldID = walberla::field::addToStorage<ScalarField_T>(concurrentBlocks, "concurrent eddy viscosity field", real_t(0), layout, fieldGhostLayers);
        BlockDataID concurrentStrainRateFieldID = walberla::field::addToStorage<SecondOrderTensorField_T>(concurrentBlocks, "concurrent strain rate field", real_t(0), layout, fieldGhostLayers);
        BlockDataID concurrentVelocityFieldID = walberla::field::addToStorage<VectorField_T>(concurrentBlocks, "concurrent velocity field", real_t(0), layout, fieldGhostLayers);
        BlockDataID concurrentFlagFieldID = walberla::field::addFlagFieldToStorage<FlagField_T>(concurrentBlocks, "concurrent flag field", fieldGhostLayers);

        const StorageSpecification_T storageSpec{};

        uint_t initialWelfordCounter = 0;
        if(precursorDomainSetup.snapshotHandler_->loadSnapshot) {
            auto & ssh = precursorDomainSetup.snapshotHandler_;

            auto vectorFieldDataHandling = std::make_shared<field::FlattenedFieldHandling<VectorField_T>> (precursorBlocks, fieldGhostLayers, layout);
            const auto forceFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->forceFile).c_str());

            precursorForceFieldID = precursorBlocks->loadBlockData(forceFileName, vectorFieldDataHandling, "force field");

            auto pdfFieldDataHandling = std::make_shared< walberla::lbm_generated::internal::PdfFieldHandling< StorageSpecification_T > >( precursorBlocks, storageSpec, fieldGhostLayers, layout );

            const auto pdfFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->pdfFile).c_str());
            precursorPdfFieldID = precursorBlocks->loadBlockData(pdfFileName, pdfFieldDataHandling, "pdf field");

            if(!ssh->meanVelFile.empty()) {
                initialWelfordCounter = ssh->welfordCounter;
                const auto meanVelFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->meanVelFile).c_str());

                precursorMeanVelocityFieldWFBID    = precursorBlocks->loadBlockData(meanVelFileName, vectorFieldDataHandling, "mean velocity field WFB");
            }
        } else {
            precursorForceFieldID = walberla::field::addToStorage<VectorField_T>(precursorBlocks, "precursor force field", real_t(0), layout, fieldGhostLayers);
            precursorPdfFieldID = walberla::lbm_generated::addPdfFieldToStorage(precursorBlocks, "precursor pdf field", storageSpec, fieldGhostLayers, layout);
            precursorMeanVelocityFieldWFBID = walberla::field::addToStorage<VectorField_T>(precursorBlocks, "precursor mean velocity field WFB", real_t(0), layout, fieldGhostLayers);
         }

        concurrentForceFieldID = walberla::field::addToStorage<VectorField_T>(concurrentBlocks, "concurrent force field", real_t(0), layout, fieldGhostLayers);
        concurrentPdfFieldID = walberla::lbm_generated::addPdfFieldToStorage(concurrentBlocks, "concurrent pdf field", storageSpec, fieldGhostLayers, layout);
        concurrentMeanVelocityFieldOutputID = walberla::field::addToStorage<VectorField_T>(concurrentBlocks, "concurrent mean velocity field output", real_t(0), layout, fieldGhostLayers);
        concurrentSumOfSquaresFieldID = walberla::field::addToStorage<SecondOrderTensorField_T>(concurrentBlocks, "sum of squares field", real_t(0), layout, fieldGhostLayers);
        concurrentMeanVelocityFieldWFBID = walberla::field::addToStorage<VectorField_T>(precursorBlocks, "precursor mean velocity field WFB", real_t(0), layout, fieldGhostLayers);

        precursorFieldMap[output::Fields::PDF] = precursorPdfFieldID;
        precursorFieldMap[output::Fields::FORCE] = precursorForceFieldID;
        precursorFieldMap[output::Fields::DENSITY] = precursorDensityFieldID;
        precursorFieldMap[output::Fields::VELOCITY] = precursorVelocityFieldID;
        precursorFieldMap[output::Fields::FLAG] = precursorFlagFieldID;

        concurrentFieldMap[output::Fields::PDF] = concurrentPdfFieldID;
        concurrentFieldMap[output::Fields::FORCE] = concurrentForceFieldID;
        concurrentFieldMap[output::Fields::DENSITY] = concurrentDensityFieldID;
        concurrentFieldMap[output::Fields::VELOCITY] = concurrentVelocityFieldID;
        concurrentFieldMap[output::Fields::MEAN_VELOCITY_OUTPUT] = concurrentMeanVelocityFieldOutputID;
        concurrentFieldMap[output::Fields::STRAIN_RATE] = concurrentStrainRateFieldID;
        concurrentFieldMap[output::Fields::EDDY_VISCOSITY] = concurrentNutFieldID;
        concurrentFieldMap[output::Fields::SUM_OF_SQUARES] = concurrentSumOfSquaresFieldID;
        concurrentFieldMap[output::Fields::FLAG] = concurrentFlagFieldID;

        /// ABL calculations
        auto initParams = globalConfig->getOneBlock("Initialisation");
        const auto initType = domain::DomainInitialisation::toType(initParams.getParameter<std::string>("type", "uniform"));
        const Vector3<real_t> initialVelocity = initParams.getParameter<Vector3<real_t>>("initialVelocity", Vector3<real_t>(0.0));

        real_t u_tau{};
        const real_t kappa = parameters.getParameter<real_t>("kappa", real_t(0.42));
        const real_t B = parameters.getParameter<real_t>("B", real_t(5.5));
        const real_t roughnessLengthRatio = parameters.getParameter<real_t>("roughnessLengthRatio", real_t(1e-4));
        const real_t referenceHeight = parameters.getParameter<real_t>("referenceHeight_LU", real_t(-1));
        const uint32_t samplingHeight = parameters.getParameter<uint32_t>("samplingHeight_LU", uint32_t(0));
        const real_t roughnessLength = roughnessLengthRatio * referenceHeight;

        const uint_t welfordIntervalOutput = walberlaEnv.config()->getOneBlock("ConcurrentOutput").getParameter<uint_t>("welfordIntervalOutput", uint_t(0));

        if(precursorBoundarySetup.inflowType() == InflowSetup::Periodic ||
            initType != domain::DomainInitialisation::UNIFORM) {
            // u_tau from MOST
            u_tau = kappa * initialVelocity[0] / std::log(real_t(1) / roughnessLengthRatio);
        }

        /// DOMAIN INITIALISATION
        WALBERLA_LOG_INFO_ON_ROOT("Initialise precursor and concurrent domain...")
        if(!precursorDomainSetup.snapshotHandler_->loadSnapshot) {

            std::unique_ptr<domain::DomainInitialisation> precursorInitialiser{};
            std::unique_ptr<domain::DomainInitialisation> concurrentInitialiser{};
            if (initType == domain::DomainInitialisation::UNIFORM) {
                WALBERLA_LOG_INFO_ON_ROOT("... Uniform")
                precursorInitialiser = std::make_unique<domain::UniformInitialisation>(precursorBlocks, initialVelocity);
                concurrentInitialiser = std::make_unique<domain::UniformInitialisation>(concurrentBlocks, initialVelocity);

            } else if (initType == domain::DomainInitialisation::ASMUTH) {
                precursorInitialiser = std::make_unique<domain::AsmuthInitialisation>(precursorBlocks, u_tau, roughnessLength, kappa, B,
                                                                             viscosity_LU, precursorDomainSetup.domainSize_, true);
                concurrentInitialiser = std::make_unique<domain::AsmuthInitialisation>(concurrentBlocks, u_tau, roughnessLength, kappa, B,
                                                                             viscosity_LU, concurrentDomainSetup.domainSize_, true);
            } else if (initType == domain::DomainInitialisation::LOG_LAW) {
                precursorInitialiser = std::make_unique<domain::AsmuthInitialisation>(precursorBlocks, u_tau, roughnessLength, kappa, B,
                                                                             viscosity_LU, precursorDomainSetup.domainSize_, false);
                concurrentInitialiser = std::make_unique<domain::AsmuthInitialisation>(concurrentBlocks, u_tau, roughnessLength, kappa, B,
                                                                             viscosity_LU, concurrentDomainSetup.domainSize_, false);
            }

            WALBERLA_LOG_INFO_ON_ROOT("... Macro setters")
            using MacroSetter_T = walberla::pystencils::waLBerlaWind_MacroSetter;

            MacroSetter_T precursorSetter(precursorDensityFieldID, precursorForceFieldID, precursorPdfFieldID, precursorVelocityFieldID);
            MacroSetter_T concurrentSetter(concurrentDensityFieldID, concurrentForceFieldID, concurrentPdfFieldID, concurrentVelocityFieldID);

            precursorInitialiser->setViaVelocityField<VectorField_T, MacroSetter_T>(precursorVelocityFieldID, precursorSetter);
            // Probably useless to initialize since we will copy it later
            concurrentInitialiser->setViaVelocityField<VectorField_T, MacroSetter_T>(concurrentVelocityFieldID, concurrentSetter);
            // WARNING: not sure if we should initialize like this here
            concurrentInitialiser->setViaVelocityField<VectorField_T, MacroSetter_T>(concurrentMeanVelocityFieldOutputID, concurrentSetter);
            concurrentInitialiser->setViaVelocityField<VectorField_T, MacroSetter_T>(concurrentMeanVelocityFieldWFBID, concurrentSetter);
            precursorInitialiser->setViaVelocityField<VectorField_T, MacroSetter_T>(precursorMeanVelocityFieldWFBID, precursorSetter);
        } else {
            if(precursorDomainSetup.snapshotHandler_->meanVelFile.empty()) {
                WALBERLA_ABORT("Mean velocity not loaded from snapshot but also not initialised.")
            }
            // initialise velocity field from pfd field if snapshot
            walberla::pystencils::waLBerlaWind_MacroGetter getter(precursorDensityFieldID, precursorForceFieldID, precursorPdfFieldID, precursorVelocityFieldID);
            for( auto & block : *precursorBlocks ) {
                getter(&block);
            }
        }

        SweepCollection_T precursorSweepCollection ( precursorBlocks, precursorDensityFieldID, precursorForceFieldID, precursorDummyNutFieldID, precursorPdfFieldID, precursorVelocityFieldID, omega );
        SweepCollection_T concurrentSweepCollection( concurrentBlocks, concurrentDensityFieldID, concurrentForceFieldID, concurrentNutFieldID, concurrentPdfFieldID, concurrentVelocityFieldID, omega );

        WALBERLA_MPI_BARRIER()
        WALBERLA_LOG_INFO_ON_ROOT("Initialisation done")

        /// BOUNDARIES

        WALBERLA_LOG_INFO_ON_ROOT("Setting up boundaries...")

        concurrentBoundarySetup.fillFlagFieldFromConfig<FlagField_T>(concurrentBlocks, concurrentFlagFieldID, FluidFlagUID,
                                                           NoSlipFlagUID, WFBFlagUID, SymmetryFlagUID, UniformInflowFlagUID, LogLawInflowFlagUID, OutflowFlagUID);
        precursorBoundarySetup.fillFlagFieldFromConfig<FlagField_T>(precursorBlocks, precursorFlagFieldID, FluidFlagUID,
                                                            NoSlipFlagUID, WFBFlagUID, SymmetryFlagUID, UniformInflowFlagUID, LogLawInflowFlagUID, OutflowFlagUID);

        if(concurrentBoundarySetup.wallType() == WallSetup::WFB){
            WALBERLA_CHECK(referenceHeight > real_t(0), "Reference height must be given in the parameter file for wall function boundary conditions.")
            WALBERLA_CHECK(concurrentBoundarySetup.environmentSetup() == EnvironmentSetup::Open,
                           "Wall-function bounce is currently only supported for an open environment")
            WALBERLA_CHECK(samplingHeight > uint_t(0), "Sampling height must be given in the parameter file for wall function boundary conditions.")
        }
        if(precursorBoundarySetup.wallType() == WallSetup::WFB){
            WALBERLA_CHECK(referenceHeight > real_t(0), "Reference height must be given in the parameter file for wall function boundary conditions.")
            WALBERLA_CHECK(precursorBoundarySetup.environmentSetup() == EnvironmentSetup::Open,
                           "Wall-function bounce is currently only supported for an open environment")
            WALBERLA_CHECK(samplingHeight > uint_t(0), "Sampling height must be given in the parameter file for wall function boundary conditions.")
        }

        const auto & u_in = concurrentBoundarySetup.inflowVelocity();
        std::function< walberla::math::Vector3< real_t >(const walberla::Cell&, const std::shared_ptr<walberla::StructuredBlockForest >&, walberla::IBlock&)>
            velocityInitialisation = std::bind(turbine_core::boundary::VelocityCallback, std::placeholders::_1,
                std::placeholders::_2, std::placeholders::_3, roughnessLength, kappa, u_tau);

        BoundaryCollection_T concurrentBoundaryCollection( concurrentBlocks, concurrentFlagFieldID, concurrentPdfFieldID, FluidFlagUID, concurrentForceFieldID,
                                                 concurrentMeanVelocityFieldWFBID, samplingHeight, roughnessLengthRatio * referenceHeight,
                                                 u_in[0], u_in[1], u_in[2], velocityInitialisation);
        BoundaryCollection_T precursorBoundaryCollection( precursorBlocks, precursorFlagFieldID, precursorPdfFieldID, FluidFlagUID, precursorForceFieldID,
                                                 precursorMeanVelocityFieldWFBID, samplingHeight, roughnessLengthRatio * referenceHeight,
                                                 u_in[0], u_in[1], u_in[2], velocityInitialisation);


        const int periodicShiftValue = precursorBoundariesConfig.getParameter<int>("periodicShiftValue", int(real_t(0.33) * real_t(precursorDomainSetup.domainSize_[2])));

        walberla::boundary::ShiftedPeriodicity<PdfField_T> shiftedPeriodicity(precursorBlocks, precursorPdfFieldID,
            fieldGhostLayers, 0, 1, periodicShiftValue);

        /// Compute inflow cells indices
        std::map<const walberla::IBlockID*, std::vector<uint_t>> cpuInflowCellsIndices;

         for (auto concurrentBlock = concurrentBlocks->begin(); concurrentBlock != concurrentBlocks->end();
              ++concurrentBlock) {

             auto concurrentFlagField = concurrentBlock->getData< FlagField_T >(concurrentFlagFieldID);
             auto concurrentPdfFieldCpu = concurrentBlock->getData< PdfField_T >(concurrentPdfFieldID);

             for (auto concurrentCellIt = concurrentFlagField->beginWithGhostLayer(); concurrentCellIt != concurrentFlagField->end(); ++concurrentCellIt)
             {
                 const walberla::Cell localConcurrentCell  = concurrentCellIt.cell();
                 if (concurrentFlagField->isPartOfMaskSet(localConcurrentCell,
                                                             concurrentFlagField->getFlag(UniformInflowFlagUID))) // ExchangeBoundaryFlagUID
                 {
                     auto inflowID =
                           (localConcurrentCell.x()) * concurrentPdfFieldCpu->xStride()
                         + (localConcurrentCell.y()) * concurrentPdfFieldCpu->yStride()
                         + (localConcurrentCell.z()) * concurrentPdfFieldCpu->zStride();
                     auto offset =
                       concurrentPdfFieldCpu->xOff()*concurrentPdfFieldCpu->xStride()
                     + concurrentPdfFieldCpu->yOff()*concurrentPdfFieldCpu->yStride()
                     + concurrentPdfFieldCpu->zOff()*concurrentPdfFieldCpu->zStride();

                     cpuInflowCellsIndices[&concurrentBlock->getId()].emplace_back(inflowID+offset);
                 }
             }
        }

        /// set IDs for farm
        farm.setBlockForest(concurrentBlocks);
        farm.setFieldIDs(concurrentDensityFieldID, concurrentVelocityFieldID, concurrentForceFieldID);


        /// WELFORD
        walberla::pystencils::waLBerlaWind_Welford welfordPrecursorSweep(precursorMeanVelocityFieldWFBID, precursorVelocityFieldID, 0.0);
        walberla::pystencils::waLBerlaWind_WelfordSOP welfordConcurrentSweepOutput(concurrentMeanVelocityFieldOutputID, concurrentSumOfSquaresFieldID, concurrentVelocityFieldID, 0.0);
        walberla::pystencils::waLBerlaWind_SoSResetter welfordConcurrentSoSResetter(concurrentSumOfSquaresFieldID);

        /// TIMELOOP

        WALBERLA_LOG_INFO_ON_ROOT("Creating time loop...")

        // create time loop
        walberla::SweepTimeloop concurrentTimeLoop( concurrentBlocks->getBlockStorage(), timesteps );
        walberla::SweepTimeloop precursorTimeLoop( precursorBlocks->getBlockStorage(), timesteps );

        farm.setTimeloop(&concurrentTimeLoop);

        /// flow driver
        walberla::wind::FlowDriverCollection precursorFlowDriver(precursorBlocks, &precursorTimeLoop, globalConfig, precursorDomainSetup, precursorForceFieldID, precursorVelocityFieldID, fieldGhostLayers);
        for( auto & block : *precursorBlocks) {
            precursorFlowDriver(&block);
        }
        walberla::wind::FlowDriverCollection concurrentFlowDriver(concurrentBlocks, &concurrentTimeLoop, globalConfig, concurrentDomainSetup, concurrentForceFieldID, concurrentVelocityFieldID, fieldGhostLayers);
        for( auto & block : *concurrentBlocks) {
            concurrentFlowDriver(&block);
        }

        // create communication for PdfField

        // create communication for fields
        walberla::blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > precursorCommunication( precursorBlocks );
        precursorCommunication.addPackInfo( std::make_shared< walberla::lbm_generated::UniformGeneratedPdfPackInfo<PdfField_T> >( precursorPdfFieldID ) );
        precursorCommunication.addPackInfo( std::make_shared< walberla::field::communication::PackInfo<ScalarField_T> >( precursorDensityFieldID ) );
        precursorCommunication.addPackInfo( std::make_shared< walberla::field::communication::PackInfo<VectorField_T > >( precursorVelocityFieldID ) );

        walberla::blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > concurrentCommunication( concurrentBlocks );
        concurrentCommunication.addPackInfo( std::make_shared< walberla::lbm_generated::UniformGeneratedPdfPackInfo<PdfField_T> >( concurrentPdfFieldID ) );
        concurrentCommunication.addPackInfo( std::make_shared< walberla::field::communication::PackInfo<ScalarField_T> >( concurrentDensityFieldID ) );
        concurrentCommunication.addPackInfo( std::make_shared< walberla::field::communication::PackInfo<VectorField_T > >( concurrentVelocityFieldID ) );

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
        precursorCommunication();
        concurrentCommunication();

        // output
        TurbineOutput_T concurrentTurbineOutput{
                concurrentBlocks, &concurrentTimeLoop,
                globalConfig->getOneBlock("ConcurrentOutput"), FluidFlagUID,
                concurrentFieldMap
        };
        TurbineOutput_T precursorTurbineOutput{
            precursorBlocks, &precursorTimeLoop,
            globalConfig->getOneBlock("PrecursorOutput"), FluidFlagUID,
            precursorFieldMap
        };

        //FIXME the current implementation expects the omega field, not the nu_t field. Strain rate writer must be called BEFORE eddy viscosity writer

        walberla::pystencils::waLBerlaWind_StrainRateWriter strainRateWriter(concurrentForceFieldID, concurrentNutFieldID, concurrentPdfFieldID, concurrentStrainRateFieldID);
        concurrentTurbineOutput.addBeforeFunction( [&]() {
            for( auto & block : *concurrentBlocks ) {
                strainRateWriter(&block);
            }});

        walberla::pystencils::waLBerlaWind_EddyViscosityWriter eddyViscosityWriter(concurrentNutFieldID, omega);
        concurrentTurbineOutput.addBeforeFunction( [&]() {
            for( auto & block : *concurrentBlocks ) {
                eddyViscosityWriter(&block);
            }});

        auto concurrentOutput = std::bind( &TurbineOutput_T::write, &concurrentTurbineOutput );
        auto precursorOutput = std::bind( &TurbineOutput_T::write, &precursorTurbineOutput );

        //NOTE must convert sweeps that are altered to lambdas, otherwise copy and counter will stay 0
        welfordPrecursorSweep.setCounter(real_t(initialWelfordCounter));
        welfordConcurrentSweepOutput.setCounter(real_t(initialWelfordCounter));
        auto welfordLambda = [&welfordPrecursorSweep](walberla::IBlock * block) {
            welfordPrecursorSweep(block);
        };
        auto welfordLambdaOutput = [&welfordConcurrentSweepOutput](walberla::IBlock * block) {
            welfordConcurrentSweepOutput(block);
        };

        concurrentTimeLoop.add() << walberla::BeforeFunction( concurrentCommunication, "Concurrent field communication" )
                                 << walberla::Sweep( concurrentBoundaryCollection.getSweep(), "Concurrent boundary handling" )
                                 << walberla::AfterFunction( [&]() {
                                     auto precursorBlock = precursorBlocks->begin();
                                     for (auto concurrentBlock = concurrentBlocks->begin(); concurrentBlock != concurrentBlocks->end();
                                         ++concurrentBlock, ++precursorBlock) {
                                         auto concurrentPDFField = concurrentBlock->getData< PdfField_T >(concurrentPdfFieldID);
                                         auto precursorPDFField  = precursorBlock->getData< PdfField_T >(precursorPdfFieldID);

                                         auto nInflowCells = cpuInflowCellsIndices[&concurrentBlock->getId()].size();
                                         copyBoundaryData( (real_t *) precursorPDFField->data(), (real_t *) concurrentPDFField->data(), cpuInflowCellsIndices[&concurrentBlock->getId()], nInflowCells, Stencil_T::Q, uint_t(concurrentPDFField->fStride()) );
                                     }
                                 }, "Concurrent boundary copy" );

        concurrentTimeLoop.add() << walberla::Sweep( concurrentSweepCollection.streamCollide() , "Concurrent LBM streamCollide sweep")
                                 << walberla::AfterFunction( concurrentOutput, "Concurrent VTK output" );

        concurrentTimeLoop.add() << walberla::BeforeFunction([&](){
                                     if(welfordIntervalOutput && (uint_t(welfordConcurrentSweepOutput.getCounter()) % welfordIntervalOutput == 0)){
                                         welfordConcurrentSweepOutput.setCounter(real_t(1));
                                         // reset the field: swap pointers
                                         // We live a dangerous life: we swap the velocity field counter
                                         for( auto blockIt = concurrentBlocks->begin(); blockIt != concurrentBlocks->end(); ++blockIt ) {
                                             auto velocityField = blockIt->getData<VectorField_T>( concurrentVelocityFieldID );
                                             auto meanVelocityField = blockIt->getData<VectorField_T>( concurrentMeanVelocityFieldOutputID );
                                             velocityField->swapDataPointers(meanVelocityField);

                                             welfordConcurrentSoSResetter(blockIt.get());
                                         }
                                     } else {
                                         welfordConcurrentSweepOutput.setCounter(real_t(welfordConcurrentSweepOutput.getCounter()+1));
                                     }
                                 }, "Concurrent welford sweep Output")
                                 << walberla::Sweep(welfordLambdaOutput, "Concurrent welford sweep Output");

        concurrentTimeLoop.add() << walberla::BeforeFunction( [&]() {farm.callback(Component::Function::UPDATE_DISCRETISATION);}, "Rotation" )
                                 << walberla::Sweep( evaluateDensityAndVelocity, "Evaluate density and velocity", sets::turbineSelector() )
                                 << walberla::AfterFunction(syncMacros, "Communicate macroscopic turbine data")
                                 << walberla::AfterFunction( applyControl , "Apply Control");

        concurrentTimeLoop.add() << walberla::Sweep(concurrentFlowDriver, "Concurrent setting driving force" );

        concurrentTimeLoop.add() << walberla::BeforeFunction( [&]() {farm.callback( Component::Function::CALCULATE_FORCES );}, "Calculate forces" )
                                 << walberla::BeforeFunction( syncForces, "Communicate force data" )
                                 << walberla::Sweep(spreadForces, "Spread forces", sets::turbineSelector() )
                                 << walberla::AfterFunction( [&]() {farm.callback(Component::Function::RESET_COMMUNICATION_DATA);}, "ResetCommunicationData" );

        precursorTimeLoop.add() << walberla::BeforeFunction( precursorCommunication, "Precursor field communication" )
                                << walberla::Sweep( precursorBoundaryCollection.getSweep(), "Precursor boundary handling" )
                                << walberla::AfterFunction( [&]() {
                                    if( precursorBoundarySetup.inflowType() == InflowSetup::ShiftedPeriodic )
                                        shiftedPeriodicity();
                                    } );
        precursorTimeLoop.add() << walberla::Sweep( precursorSweepCollection.stream() , "Precursor LBM stream sweep");
        precursorTimeLoop.add() << walberla::Sweep( precursorSweepCollection.collide() , "Precursor LBM collide sweep")
                                << walberla::AfterFunction( precursorOutput, "Precursor VTK output" )
                                << walberla::AfterFunction( [&](){
                                    precursorDomainSetup.snapshotHandler_->createSnapshot(precursorTimeLoop.getCurrentTimeStep(),
                                        precursorBlocks, precursorPdfFieldID, precursorForceFieldID, uint_t(welfordPrecursorSweep.getCounter()),
                                        precursorMeanVelocityFieldWFBID);
                                }, "Snapshot creation");

        precursorTimeLoop.add() << walberla::BeforeFunction([&](){welfordPrecursorSweep.setCounter(real_t(welfordPrecursorSweep.getCounter()+1));}, "Precursor welford sweep")
                                << walberla::Sweep(welfordLambda, "welford sweep");
        precursorTimeLoop.add() << walberla::Sweep(precursorFlowDriver, "Precursor setting driving force" );

        // LBM stability check
        concurrentTimeLoop.addFuncAfterTimeStep( walberla::makeSharedFunctor( walberla::field::makeStabilityChecker< PdfField_T, FlagField_T >(globalConfig, concurrentBlocks, concurrentPdfFieldID, concurrentFlagFieldID, FluidFlagUID ) ),
                                       "Concurrent LBM stability check" );

        // log remaining time
        precursorTimeLoop.addFuncAfterTimeStep( walberla::timing::RemainingTimeLogger( precursorTimeLoop.getNrOfTimeSteps(), remainingTimeLoggerFrequency ), "remaining time logger" );

        WALBERLA_LOG_INFO_ON_ROOT("Running timeloop...")

        WALBERLA_MPI_WORLD_BARRIER()

        walberla::WcTimingPool concurrentTiming;
        walberla::WcTimingPool precursorTiming;
        concurrentTiming.registerTimer("Force Output");

        walberla::WcTimer concurrentTimer;
        walberla::WcTimer precursorTimer;
        precursorTimer.start();
        for( uint_t i = 0; i < timesteps; ++i ) {
            precursorTimeLoop.singleStep(precursorTiming);

            // Copy pdf fields
            if (i == concurrentStart)
            {
                WALBERLA_LOG_INFO_ON_ROOT("Concurrent starts")

                concurrentTimer.start();

                auto concurrentBlock = concurrentBlocks->begin();
                for (auto precursorBlock = precursorBlocks->begin(); precursorBlock != precursorBlocks->end();
                     ++precursorBlock)
                {
                    const auto gl = cell_idx_t(fieldGhostLayers);
                    auto precursorPDFField  = precursorBlock->getData< PdfField_T >(precursorPdfFieldID);
                    auto concurrentPDFField   = concurrentBlock->getData< PdfField_T >(concurrentPdfFieldID);
                    std::copy( &precursorPDFField->get(-gl,-gl,-gl, 0),
                        &precursorPDFField->get(+gl,+gl,+gl, PdfField_T::F_SIZE-1),
                        &concurrentPDFField->get(-gl,-gl,-gl, 0) );

                    auto precursorDensityField  = precursorBlock->getData< ScalarField_T >(precursorDensityFieldID);
                    auto concurrentDensityField   = concurrentBlock->getData< ScalarField_T >(concurrentDensityFieldID);
                    std::copy( &precursorDensityField->get(-gl,-gl,-gl, 0),
                        &precursorDensityField->get(+gl,+gl,+gl, ScalarField_T::F_SIZE-1),
                        &concurrentDensityField->get(-gl,-gl,-gl, 0) );

                    auto precursorVelocityField  = precursorBlock->getData< VectorField_T >(precursorVelocityFieldID);
                    auto concurrentVelocityField   = concurrentBlock->getData< VectorField_T >(concurrentVelocityFieldID);
                    std::copy( &precursorVelocityField->get(-gl,-gl,-gl, 0),
                        &precursorVelocityField->get(+gl,+gl,+gl, VectorField_T::F_SIZE-1),
                        &concurrentVelocityField->get(-gl,-gl,-gl, 0) );

                    auto precursorMeanVelocityField  = precursorBlock->getData< VectorField_T >(precursorMeanVelocityFieldWFBID);
                    auto concurrentMeanVelocityField   = concurrentBlock->getData< VectorField_T >(concurrentMeanVelocityFieldWFBID);
                    std::copy( &precursorMeanVelocityField->get(-gl,-gl,-gl, 0),
                        &precursorMeanVelocityField->get(+gl,+gl,+gl, VectorField_T::F_SIZE-1),
                        &concurrentMeanVelocityField->get(-gl,-gl,-gl, 0) );

                    concurrentBlock++;
                }
            }

            if (i >= concurrentStart)
            {
                // Copy the mean velocity field from precursor to concurrent
                auto concurrentBlock = concurrentBlocks->begin();
                for (auto precursorBlock = precursorBlocks->begin(); precursorBlock != precursorBlocks->end();
                     ++precursorBlock) {
                    auto precursorMeanVelocityFieldWFB  = precursorBlock->getData< VectorField_T >(precursorMeanVelocityFieldWFBID);
                    auto concurrentMeanVelocityFieldWFB  = concurrentBlock->getData< VectorField_T >(concurrentMeanVelocityFieldWFBID);
                    concurrentMeanVelocityFieldWFB->set(precursorMeanVelocityFieldWFB);
                    concurrentBlock++;
                }

                concurrentTimeLoop.singleStep(concurrentTiming);
            }

            concurrentTiming["Force Output"].start();
            // farm.callback(Component::Output::ORIENTATIONS);
            for(auto & block : *concurrentBlocks) {
                farm.writeForceOutput(&block, forceOutputStart, forceOutputFrequency);
            }
            farm.writeTurbinePowerAndThrust();
            concurrentTiming["Force Output"].end();
        }

        concurrentTimer.end();
        precursorTimer.end();

        double concurrentTime = concurrentTimer.max();
        walberla::mpi::reduceInplace(concurrentTime, walberla::mpi::MAX);
        double precursorTime = precursorTimer.max();
        walberla::mpi::reduceInplace(precursorTime, walberla::mpi::MAX);

        const auto concurrentTimeloopTiming = concurrentTiming.getReduced();
        WALBERLA_LOG_INFO_ON_ROOT("Concurrent timeloop timing:\n" << *concurrentTimeloopTiming)
        const auto precursorTimeloopTiming = precursorTiming.getReduced();
        WALBERLA_LOG_INFO_ON_ROOT("Precursor timeloop timing:\n" << *precursorTimeloopTiming)

        walberla::lbm::PerformanceEvaluation<FlagField_T> concurrentPerformance(concurrentBlocks, concurrentFlagFieldID, FluidFlagUID);
        concurrentPerformance.logResultOnRoot(timesteps, concurrentTime);

        walberla::lbm::PerformanceEvaluation<FlagField_T> precursorPerformance(precursorBlocks, precursorFlagFieldID, FluidFlagUID);
        precursorPerformance.logResultOnRoot(timesteps, precursorTime);

        farm.calculatePowerAndThrust();

        logging::WindTurbineApplicationLogger<domain::DomainSetup, domain::BoundarySetup, codegen::KernelInfo, walberla::timing::WcPolicy> concurrentLogger {
            concurrentDomainSetup, concurrentBoundarySetup, windFarmConfig, &concurrentTiming, concurrentPerformance.loggingString(timesteps, concurrentTime),
            false, "WindTurbine_Logging_concurrent.log"
        };
        logging::WindTurbineApplicationLogger<domain::DomainSetup, domain::BoundarySetup, codegen::KernelInfo, walberla::timing::WcPolicy> precursorLogger {
            precursorDomainSetup, precursorBoundarySetup, windFarmConfig, &precursorTiming, precursorPerformance.loggingString(timesteps, precursorTime),
            false, "WindTurbine_Logging_precursor.log"
        };

        return EXIT_SUCCESS;

    }

}

int main(int argc, char** argv) {
    return turbine_core::main(argc, argv);
}
