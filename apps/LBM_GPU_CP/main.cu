#include <iostream>
#include <cmath>
#include <vector>
#include <functional>

#include <core/Environment.h>

#include "output/Logging.h"

//TODO
#include <blockforest/Initialization.h>
#include <core/all.h>
#include <core/mpi/all.h>
#include <domain_decomposition/all.h>
#include <field/all.h>
#include <lbm_generated/field/AddToStorage.h>
#include <lbm_generated/gpu/UniformGeneratedGPUPdfPackInfo.h>
#include <lbm/all.h>
#include <lbm/vtk/all.h>
#include <timeloop/all.h>
#include <gpu/AddGPUFieldToStorage.h>
#include <gpu/HostFieldAllocator.h>
#include <gpu/communication/UniformGPUScheme.h>
#include <gpu/communication/MemcpyPackInfo.h>
#include <gpu/ShiftedPeriodicity.h>

#include "walberla_helper/field/all.h"
#include "walberla_helper/field/FieldHandling.h"

#include "domain/BoundarySetup.h"
#include "domain/BoundaryHandling.h"
#include "domain/DomainSetup.h"
#include "domain/DomainInitialiser.h"

#include "external_forces/FlowDriverGPU.h"

#include "farm_creator/WindFarm.h"
#include "farm_creator/ConfigTurbineCreator.h"

#include "wind_turbine_core/ProjectDefines.h"

#include "conversion/Conversion.h"

#include "turbine_topology/gpu/Tree.h"

#include "WindTurbineDefines.h"

#include "output/TurbineOutput.h"

#include "waLBerlaWind_KernelInfo.h"

namespace turbine_core {

    using InterpolationKernel_T = projectors::SmoothedDiracDeltaKernel;
    using DistributionKernel_T = projectors::SmoothedDiracDeltaKernel;

    using WindFarm_T = WindFarm<topology::gpu::Tree, creator::ConfigTurbineCreator,
            ScalarField_T, VectorField_T, VectorField_T,
            //            TrilinearMacroscopicFieldInterpolator,
            KernelFieldInterpolator_T<InterpolationKernel_T>,
            KernelDistributor_T<DistributionKernel_T>>;

    using TurbineOutput_T = output::TurbineOutput<
            FlagField_T, ScalarField_T, VectorField_T, SecondOrderTensorField_T, ThirdOrderTensorField_T,
            PdfField_T, Stencil_T, StorageSpecification_T::zeroCenteredPDFs, StorageSpecification_T::compressible
    >;

    template<typename Type_T>
    using GPUField_T = walberla::gpu::GPUField<Type_T>;

    GLOBAL_PREFIX void copyBoundaryData( real_t * precursorField, real_t * concurrentField,
        uint_t * cellIndices, const uint_t nbCells, const uint_t q, const uint_t fStride ) {
        // Make sure the field is in fxyz layout
        uint_t id = blockIdx.x*blockDim.x+threadIdx.x;

        if( id < nbCells ) {
            for (uint_t qIdx = 0; qIdx < q; ++qIdx) {
                concurrentField[cellIndices[id] + qIdx*fStride] = precursorField[cellIndices[id] + qIdx*fStride];
            }
        }
    }

    GLOBAL_PREFIX void setFieldFZYX( real_t * data, const int fStride,
                                     const real_t xValue, const real_t yValue, const real_t zValue ) {

        uint_t id = blockIdx.x*blockDim.x+threadIdx.x;

        if( id < fStride ) {
            data[id              ] = xValue;
            data[id +     fStride] = yValue;
            data[id + 2 * fStride] = zValue;
        }

    }

    int main(int argc, char** argv) {

        walberla::Environment walberlaEnv(argc, argv);

        walberla::logging::Logging::instance()->setLogLevel(walberla::logging::Logging::INFO);

        /// general simulation parameters
        auto globalConfig = walberlaEnv.config();
        auto parameters = globalConfig->getOneBlock("Parameters");

        uint_t timesteps = parameters.getParameter<uint_t>("timesteps", uint_t(10)); ++timesteps;
        const uint_t concurrentStart = parameters.getParameter<uint_t>("concurrentStart", uint_t(5000));
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

        auto precursorBoundariesConfig = globalConfig->getOneBlock( "PrecursorBoundaries" );
        domain::BoundarySetup precursorBoundarySetup(precursorBoundariesConfig);
        auto concurrentBoundariesConfig = globalConfig->getOneBlock( "ConcurrentBoundaries" );
        domain::BoundarySetup concurrentBoundarySetup(concurrentBoundariesConfig);

        domain::DomainSetup precursorDomainSetup ( globalConfig, precursorBoundarySetup.periodicity() );
        auto precursorBlocks = precursorDomainSetup.createUniformBlockForest(farm.getTurbineAABBs());
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

        //CPU
        BlockDataID precursorPdfFieldCpuID;
        BlockDataID precursorMeanVelocityFieldWFBCpuID = walberla::field::addToStorage<VectorField_T>(precursorBlocks, "Precursor mean WFB velocity field CPU", real_t(0), layout, fieldGhostLayers);

        BlockDataID precursorDensityFieldCpuID = walberla::field::addToStorage<ScalarField_T>(precursorBlocks, "Precursor density field CPU", real_t(1), layout, fieldGhostLayers);
        BlockDataID precursorVelocityFieldCpuID = walberla::field::addToStorage<VectorField_T>(precursorBlocks, "Precursor velocity field CPU", real_t(0), layout, fieldGhostLayers);
        BlockDataID precursorForceFieldCpuID    = walberla::field::addToStorage<VectorField_T>(precursorBlocks, "Precursor force field CPU", real_t(0), layout, fieldGhostLayers);
        BlockDataID precursorDummyNutFieldCpuID = walberla::field::addToStorage<ScalarField_T>(precursorBlocks, "Precursor eddy viscosity field CPU", real_t(0), layout, fieldGhostLayers);
        BlockDataID precursorFlagFieldCpuID = walberla::field::addFlagFieldToStorage<FlagField_T>(precursorBlocks, "Precursor flag field", fieldGhostLayers);
        BlockDataID precursorDummySumOfSquaresFieldCpuID = walberla::field::addToStorage<SecondOrderTensorField_T>(precursorBlocks, "Precursor sum of squares field CPU", real_t(0), layout, fieldGhostLayers);

        BlockDataID concurrentPdfFieldCpuID;
        BlockDataID concurrentMeanVelocityFieldOutputCpuID;
        BlockDataID concurrentMeanVelocityFieldWFBCpuID = walberla::field::addToStorage<VectorField_T>(concurrentBlocks, "Concurrent mean WFB velocity field CPU", real_t(0), layout, fieldGhostLayers);
        BlockDataID concurrentDensityFieldCpuID = walberla::field::addToStorage<ScalarField_T>(concurrentBlocks, "Concurrent density field CPU", real_t(1), layout, fieldGhostLayers);
        BlockDataID concurrentVelocityFieldCpuID = walberla::field::addToStorage<VectorField_T>(concurrentBlocks, "Concurrent velocity field CPU", real_t(0), layout, fieldGhostLayers);
        BlockDataID concurrentForceFieldCpuID    = walberla::field::addToStorage<VectorField_T>(concurrentBlocks, "Concurrent force field CPU", real_t(0), layout, fieldGhostLayers);
        BlockDataID concurrentNutFieldCpuID = walberla::field::addToStorage<ScalarField_T>(concurrentBlocks, "Concurrent eddy viscosity field CPU", real_t(0), layout, fieldGhostLayers);
        BlockDataID concurrentFlagFieldCpuID = walberla::field::addFlagFieldToStorage<FlagField_T>(concurrentBlocks, "Concurrent flag field", fieldGhostLayers);
        BlockDataID concurrentStrainRateFieldCpuID = walberla::field::addToStorage<SecondOrderTensorField_T>(concurrentBlocks, "strain rate field CPU", real_t(0), layout, fieldGhostLayers);
        BlockDataID concurrentSumOfSquaresFieldCpuID = walberla::field::addToStorage<SecondOrderTensorField_T>(concurrentBlocks, "sum of squares field CPU", real_t(0), layout, fieldGhostLayers);


        const StorageSpecification_T storageSpec{};
        auto allocator = std::make_shared< walberla::gpu::HostFieldAllocator<real_t> >();

        // For the snapshot handler, only the precursor cases are delt with.
        // The concurrent restart could be interesting in the future, but not restart system for the wind turbines are available yet,
        // implying we must let the wake develop again.
        uint_t initialWelfordCounter = 0;
        if(precursorDomainSetup.snapshotHandler_->loadSnapshot) {
            auto & ssh = precursorDomainSetup.snapshotHandler_;

            auto vectorFieldDataHandling = std::make_shared<field::FlattenedFieldHandling<VectorField_T>> (precursorBlocks, fieldGhostLayers, layout, allocator);
            const auto forceFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->forceFile).c_str());

            precursorForceFieldCpuID = precursorBlocks->loadBlockData(forceFileName, vectorFieldDataHandling, "force field");

            auto pdfFieldDataHandling = std::make_shared< walberla::lbm_generated::internal::PdfFieldHandling< StorageSpecification_T > >( precursorBlocks, storageSpec, fieldGhostLayers, layout, allocator );

            const auto pdfFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->pdfFile).c_str());
            precursorPdfFieldCpuID = precursorBlocks->loadBlockData(pdfFileName, pdfFieldDataHandling, "pdf field");


            if(!ssh->meanVelFile.empty()) {
                initialWelfordCounter = ssh->welfordCounter;
                const auto meanVelFileName = std::string(walberla::filesystem::path(ssh->baseFolder) / walberla::filesystem::path(ssh->meanVelFile).c_str());

                precursorMeanVelocityFieldWFBCpuID = precursorBlocks->loadBlockData(meanVelFileName, vectorFieldDataHandling, "mean velocity field output");
            }
        } else {
            const auto pdfAlloc = std::make_shared<walberla::field::AllocateAligned<real_t,256>>();
            const auto empty = walberla::Set<walberla::UID<walberla::uid::suidgenerator::S>>::emptySet();

            precursorForceFieldCpuID = walberla::field::addToStorage<VectorField_T>(precursorBlocks, "precursor force field", real_t(0), layout, fieldGhostLayers);
            precursorPdfFieldCpuID = walberla::lbm_generated::addPdfFieldToStorage(precursorBlocks, "precursor pdf field", storageSpec, fieldGhostLayers, layout,
                empty, empty, pdfAlloc);
            precursorMeanVelocityFieldWFBCpuID = walberla::field::addToStorage<VectorField_T>(precursorBlocks, "precursor mean velocity field", real_t(0), layout, fieldGhostLayers);
            precursorDummySumOfSquaresFieldCpuID = walberla::field::addToStorage<SecondOrderTensorField_T>(precursorBlocks, "precursor sum of squares field", real_t(0), layout, fieldGhostLayers);

            concurrentForceFieldCpuID = walberla::field::addToStorage<VectorField_T>(concurrentBlocks, "concurrent force field", real_t(0), layout, fieldGhostLayers);
            concurrentPdfFieldCpuID = walberla::lbm_generated::addPdfFieldToStorage(concurrentBlocks, "concurrent pdf field", storageSpec, fieldGhostLayers, layout,
                empty, empty, pdfAlloc);
            concurrentMeanVelocityFieldOutputCpuID = walberla::field::addToStorage<VectorField_T>(concurrentBlocks, "concurrent mean velocity field output", real_t(0), layout, fieldGhostLayers);
            concurrentMeanVelocityFieldWFBCpuID = walberla::field::addToStorage<VectorField_T>(concurrentBlocks, "concurrent mean velocity field WFB", real_t(0), layout, fieldGhostLayers);
            concurrentSumOfSquaresFieldCpuID = walberla::field::addToStorage<SecondOrderTensorField_T>(concurrentBlocks, "sum of squares field", real_t(0), layout, fieldGhostLayers);
        }

        precursorFieldMap[output::Fields::PDF] = precursorPdfFieldCpuID;
        precursorFieldMap[output::Fields::FORCE] = precursorForceFieldCpuID;
        precursorFieldMap[output::Fields::VELOCITY] = precursorVelocityFieldCpuID;
        precursorFieldMap[output::Fields::MEAN_VELOCITY_WFB] = precursorMeanVelocityFieldWFBCpuID;
        precursorFieldMap[output::Fields::FLAG] = precursorFlagFieldCpuID;

        concurrentFieldMap[output::Fields::PDF] = concurrentPdfFieldCpuID;
        concurrentFieldMap[output::Fields::FORCE] = concurrentForceFieldCpuID;
        concurrentFieldMap[output::Fields::VELOCITY] = concurrentVelocityFieldCpuID;
        concurrentFieldMap[output::Fields::MEAN_VELOCITY_OUTPUT] = concurrentMeanVelocityFieldOutputCpuID;
        concurrentFieldMap[output::Fields::STRAIN_RATE] = concurrentStrainRateFieldCpuID;
        concurrentFieldMap[output::Fields::EDDY_VISCOSITY] = concurrentNutFieldCpuID;
        concurrentFieldMap[output::Fields::SUM_OF_SQUARES] = concurrentSumOfSquaresFieldCpuID;
        concurrentFieldMap[output::Fields::FLAG] = concurrentFlagFieldCpuID;

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

        const uint_t welfordIntervalOutput = walberlaEnv.config()->getOneBlock("ConcurrentOutput").getParameter<uint_t>("welfordIntervalOutput", uint_t(0));

        const real_t u_tau = kappa * initialVelocity[0] / log(real_t(1) / roughnessLengthRatio);

        /// DOMAIN INITIALISATION

        WALBERLA_LOG_INFO_ON_ROOT("Initialise domain...")
        if(!precursorDomainSetup.snapshotHandler_->loadSnapshot) {

            std::unique_ptr<domain::DomainInitialisation> precursorInitialiser{};
            std::unique_ptr<domain::DomainInitialisation> concurrentInitialiser{};
            if (initType == domain::DomainInitialisation::UNIFORM) {
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

            using MacroSetter_T = walberla::pystencils::waLBerlaWind_MacroSetter;
            MacroSetter_T precursorSetter(precursorDensityFieldCpuID, precursorForceFieldCpuID, precursorPdfFieldCpuID, precursorVelocityFieldCpuID);
            MacroSetter_T concurrentSetter(concurrentDensityFieldCpuID, concurrentForceFieldCpuID, concurrentPdfFieldCpuID, concurrentVelocityFieldCpuID);

            precursorInitialiser->setViaVelocityField<VectorField_T, MacroSetter_T>(precursorVelocityFieldCpuID, precursorSetter);
            precursorInitialiser->setViaVelocityField<VectorField_T, MacroSetter_T>(precursorMeanVelocityFieldWFBCpuID, precursorSetter);

            concurrentInitialiser->setViaVelocityField<VectorField_T, MacroSetter_T>(concurrentVelocityFieldCpuID, concurrentSetter);
            concurrentInitialiser->setViaVelocityField<VectorField_T, MacroSetter_T>(concurrentMeanVelocityFieldOutputCpuID, concurrentSetter);

            // setter sweep only initializes interior of domain - for push schemes to work a first communication is required here
            walberla::blockforest::communication::UniformBufferedScheme<CommunicationStencil_T> concurrentInitialComm(concurrentBlocks);
            concurrentInitialComm.addPackInfo( std::make_shared< walberla::field::communication::PackInfo<PdfField_T> >( concurrentPdfFieldCpuID ) );
            concurrentInitialComm();

            walberla::blockforest::communication::UniformBufferedScheme<CommunicationStencil_T> precursorInitialComm(precursorBlocks);
            precursorInitialComm.addPackInfo( std::make_shared< walberla::field::communication::PackInfo<PdfField_T> >( precursorPdfFieldCpuID ) );
            precursorInitialComm();

        } else {
            if(precursorDomainSetup.snapshotHandler_->meanVelFile.empty()) {
                WALBERLA_ABORT("Mean velocity not loaded from snapshot but also not initialised.")
            }
            // initialise velocity field from pfd field if snapshot
            walberla::pystencils::waLBerlaWind_MacroGetter getter(precursorDensityFieldCpuID, precursorForceFieldCpuID, precursorPdfFieldCpuID, precursorVelocityFieldCpuID);
            for( auto & block : *precursorBlocks ) {
                getter(&block);
            }
        }

        WALBERLA_LOG_INFO_ON_ROOT("Add GPU fields...")

        // GPU
        BlockDataID precursorPdfFieldGpuID = walberla::lbm_generated::addGPUPdfFieldToStorage<PdfField_T, StorageSpecification_T >(precursorBlocks, precursorPdfFieldCpuID, storageSpec, "pdf field");
        BlockDataID precursorForceFieldGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T >( precursorBlocks, precursorForceFieldCpuID, "force field GPU", true );
        BlockDataID precursorDensityFieldGpuID = walberla::gpu::addGPUFieldToStorage<ScalarField_T >( precursorBlocks, precursorDensityFieldCpuID, "density field GPU", true );

        BlockDataID precursorVelocityFieldGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T >( precursorBlocks, precursorVelocityFieldCpuID, "velocity field GPU", true );
        BlockDataID precursorMeanVelocityFieldWFBGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T >( precursorBlocks, precursorMeanVelocityFieldWFBCpuID, "meanVelocityWFB field GPU", true );
        BlockDataID precursorDummyMeanVelocityFieldOutputGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T >( precursorBlocks, precursorMeanVelocityFieldWFBCpuID, "meanVelocityOutput field GPU", true );
        BlockDataID precursorDummyNutFieldGpuID = walberla::gpu::addGPUFieldToStorage<ScalarField_T >( precursorBlocks, precursorDummyNutFieldCpuID, "nut field GPU", true );

        BlockDataID concurrentPdfFieldGpuID = walberla::lbm_generated::addGPUPdfFieldToStorage<PdfField_T, StorageSpecification_T >(concurrentBlocks, concurrentPdfFieldCpuID, storageSpec, "pdf field");
        BlockDataID concurrentDensityFieldGpuID = walberla::gpu::addGPUFieldToStorage<ScalarField_T >( concurrentBlocks, concurrentDensityFieldCpuID, "density field GPU", true );
        BlockDataID concurrentForceFieldGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T >( concurrentBlocks, concurrentForceFieldCpuID, "force field GPU", true );
        BlockDataID concurrentVelocityFieldGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T >( concurrentBlocks, concurrentVelocityFieldCpuID, "velocity field GPU", true );
        BlockDataID concurrentMeanVelocityFieldWFBGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T >( concurrentBlocks, concurrentVelocityFieldCpuID, "meanVelocityWFB field GPU", true );
        BlockDataID concurrentMeanVelocityFieldOutputGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T >( concurrentBlocks, concurrentMeanVelocityFieldOutputCpuID, "meanVelocityOutput field GPU", true );
        BlockDataID concurrentNutFieldGpuID = walberla::gpu::addGPUFieldToStorage<ScalarField_T >( concurrentBlocks, concurrentNutFieldCpuID, "nut field GPU", true );
        BlockDataID concurrentStrainRateFieldGpuID = walberla::gpu::addGPUFieldToStorage<SecondOrderTensorField_T >( concurrentBlocks, concurrentStrainRateFieldCpuID, "strainRateField field GPU", true );
        BlockDataID concurrentSumOfSquaresFieldGpuID = walberla::gpu::addGPUFieldToStorage<SecondOrderTensorField_T >(concurrentBlocks, concurrentSumOfSquaresFieldCpuID, "sumO of squares field GPU", true );

        Vector3 <int32_t> gpuBlockSize = parameters.getParameter<Vector3<int32_t>>("gpuBlockSize", Vector3<int32_t>(256, 1, 1));
        SweepCollection_T precursorSweepCollection( precursorBlocks, precursorDensityFieldGpuID, precursorForceFieldGpuID, precursorDummyNutFieldGpuID, precursorPdfFieldGpuID, precursorVelocityFieldGpuID, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2], omega );
        SweepCollection_T concurrentSweepCollection( concurrentBlocks, concurrentDensityFieldGpuID, concurrentForceFieldGpuID, concurrentNutFieldGpuID, concurrentPdfFieldGpuID, concurrentVelocityFieldGpuID, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2], omega );

        WALBERLA_MPI_BARRIER()
        WALBERLA_LOG_INFO_ON_ROOT("Initialisation done")

        /// BOUNDARIES

        WALBERLA_LOG_INFO_ON_ROOT("Setting up boundaries...")

        precursorBoundarySetup.fillFlagFieldFromConfig<FlagField_T>(precursorBlocks, precursorFlagFieldCpuID, FluidFlagUID,
                                                           NoSlipFlagUID, WFBFlagUID, SymmetryFlagUID, UniformInflowFlagUID, LogLawInflowFlagUID, OutflowFlagUID);
        concurrentBoundarySetup.fillFlagFieldFromConfig<FlagField_T>(concurrentBlocks, concurrentFlagFieldCpuID, FluidFlagUID,
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

        BoundaryCollection_T concurrentBoundaryCollection( concurrentBlocks, concurrentFlagFieldCpuID, concurrentPdfFieldGpuID, FluidFlagUID, concurrentForceFieldGpuID,
                                                 concurrentMeanVelocityFieldWFBGpuID, samplingHeight, roughnessLengthRatio * referenceHeight,
                                                 u_in[0], u_in[1], u_in[2], velocityInitialisation, concurrentPdfFieldCpuID);
        BoundaryCollection_T precursorBoundaryCollection( precursorBlocks, precursorFlagFieldCpuID, precursorPdfFieldGpuID, FluidFlagUID, precursorForceFieldGpuID,
                                                 precursorMeanVelocityFieldWFBGpuID, samplingHeight, roughnessLengthRatio * referenceHeight,
                                                 u_in[0], u_in[1], u_in[2], velocityInitialisation, precursorPdfFieldCpuID);

        const int periodicShiftValue = precursorBoundariesConfig.getParameter<int>("periodicShiftValue", int(real_t(0.33) * real_t(precursorDomainSetup.domainSize_[2])));
        walberla::gpu::ShiftedPeriodicityGPU<PdfFieldGPU_T> shiftedPeriodicity(precursorBlocks, precursorPdfFieldGpuID, fieldGhostLayers,
                                                                               0, 1, periodicShiftValue);

        /// Compute inflow cells indices
        std::map<const walberla::IBlockID*, std::vector<uint_t>> inflowCellsIndices;
        std::map<const walberla::IBlockID*, uint_t*> gpuInflowCellsIndices;

         for (auto concurrentBlock = concurrentBlocks->begin(); concurrentBlock != concurrentBlocks->end();
              ++concurrentBlock) {

             auto concurrentFlagFieldCpu = concurrentBlock->getData< FlagField_T >(concurrentFlagFieldCpuID);
             auto concurrentPdfFieldCpu = concurrentBlock->getData< PdfField_T >(concurrentPdfFieldCpuID);

             for (auto concurrentCellIt = concurrentFlagFieldCpu->beginWithGhostLayer(); concurrentCellIt != concurrentFlagFieldCpu->end(); ++concurrentCellIt)
             {
                 const walberla::Cell localConcurrentCell  = concurrentCellIt.cell();
                 if (concurrentFlagFieldCpu->isPartOfMaskSet(localConcurrentCell,
                                                             concurrentFlagFieldCpu->getFlag(UniformInflowFlagUID))) // ExchangeBoundaryFlagUID
                 {
                     auto inflowID =
                           (localConcurrentCell.x()) * concurrentPdfFieldCpu->xStride()
                         + (localConcurrentCell.y()) * concurrentPdfFieldCpu->yStride()
                         + (localConcurrentCell.z()) * concurrentPdfFieldCpu->zStride();
                     auto offset =
                       concurrentPdfFieldCpu->xOff()*concurrentPdfFieldCpu->xStride()
                     + concurrentPdfFieldCpu->yOff()*concurrentPdfFieldCpu->yStride()
                     + concurrentPdfFieldCpu->zOff()*concurrentPdfFieldCpu->zStride();

                     inflowCellsIndices[&concurrentBlock->getId()].emplace_back(inflowID+offset);
                 }
             }
             uint_t * tmpIndicesVector;
             auto blockId = &(concurrentBlock->getId());
             cudaMalloc((void**)&tmpIndicesVector, inflowCellsIndices[blockId].size() * sizeof(uint_t));
             cudaMemcpy(tmpIndicesVector, inflowCellsIndices[blockId].data(), inflowCellsIndices[blockId].size() * sizeof(uint_t), cudaMemcpyHostToDevice);
             gpuInflowCellsIndices[blockId]  = tmpIndicesVector;
             TURBINE_GPU_CHECK()
        }

        /// set IDs for farm
        farm.setBlockForest(concurrentBlocks);
        farm.setFieldIDs(precursorDensityFieldGpuID, precursorVelocityFieldGpuID, precursorForceFieldGpuID);
        farm.setFieldIDs(concurrentDensityFieldGpuID, concurrentVelocityFieldGpuID, concurrentForceFieldGpuID);

        walberla::pystencils::waLBerlaWind_SoSResetter sopResetter(concurrentSumOfSquaresFieldGpuID);

        /// WELFORD

        walberla::pystencils::waLBerlaWind_Welford precursorWelfordSweepWFB(precursorMeanVelocityFieldWFBGpuID, precursorVelocityFieldGpuID, 0.0, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);
        walberla::pystencils::waLBerlaWind_Welford concurrentWelfordSweepWFB(concurrentMeanVelocityFieldWFBGpuID, concurrentVelocityFieldGpuID, 0.0, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);
        walberla::pystencils::waLBerlaWind_WelfordSOP concurrentWelfordSweepOutput(concurrentMeanVelocityFieldOutputGpuID, concurrentSumOfSquaresFieldGpuID, concurrentVelocityFieldGpuID, 0.0, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);

        /// GPU communication

        WALBERLA_LOG_INFO_ON_ROOT("GPU communication...")

        // Communication setup
        bool cudaEnabledMPI = parameters.getParameter<bool>( "cudaEnabledMPI", false );
        if(walberla::mpi::MPIManager::instance()->numProcesses() == 1)
            cudaEnabledMPI = false;

        walberla::gpu::communication::UniformGPUScheme< CommunicationStencil_T > precursorCommunication(precursorBlocks, cudaEnabledMPI);
        precursorCommunication.addPackInfo(std::make_shared<walberla::lbm_generated::UniformGeneratedGPUPdfPackInfo< PdfFieldGPU_T >>(precursorPdfFieldGpuID));
        precursorCommunication.addPackInfo(std::make_shared<walberla::gpu::communication::MemcpyPackInfo< GPUField_T<real_t> >>(precursorDensityFieldGpuID));
        precursorCommunication.addPackInfo(std::make_shared<walberla::gpu::communication::MemcpyPackInfo< GPUField_T<real_t> >>(precursorVelocityFieldGpuID));

        walberla::gpu::communication::UniformGPUScheme< CommunicationStencil_T > concurrentCommunication(concurrentBlocks, cudaEnabledMPI);
        concurrentCommunication.addPackInfo(std::make_shared<walberla::lbm_generated::UniformGeneratedGPUPdfPackInfo< PdfFieldGPU_T >>(concurrentPdfFieldGpuID));
        concurrentCommunication.addPackInfo(std::make_shared<walberla::gpu::communication::MemcpyPackInfo< GPUField_T<real_t> >>(concurrentDensityFieldGpuID));
        concurrentCommunication.addPackInfo(std::make_shared<walberla::gpu::communication::MemcpyPackInfo< GPUField_T<real_t> >>(concurrentVelocityFieldGpuID));

        Vector3<uint_t> innerOuterSplit = parameters.getParameter<Vector3<uint_t>>("innerOuterSplit", Vector3<uint_t>(1, 1, 1));
        for(uint_t i = 0; i < 3; ++i) {
            if( int(concurrentBlocks->getNumberOfCellsPerBlock(i)) <= innerOuterSplit[i] * 2) {
                WALBERLA_ABORT_NO_DEBUG_INFO("innerOuterSplit too large - make it smaller or increase cellsPerBlock")
            }
        }

        int streamHighPriority = 0;
        int streamLowPriority = 0;
        WALBERLA_GPU_CHECK( cudaDeviceGetStreamPriorityRange(&streamLowPriority, &streamHighPriority) )
        WALBERLA_CHECK(gpuBlockSize[2] == 1)

        auto defaultStream = walberla::gpu::StreamRAII::newPriorityStream( streamLowPriority );
        auto innerOuterStreams = walberla::gpu::ParallelStreams( streamHighPriority );
        auto boundaryOuterStreams = walberla::gpu::ParallelStreams( streamHighPriority );
        auto boundaryInnerStreams = walberla::gpu::ParallelStreams( streamHighPriority );

        /// TIMELOOP

        WALBERLA_LOG_INFO_ON_ROOT("Creating time loop...")

        // create time loops
#ifdef NDEBUG
        using TimingPolicy_T = walberla::timing::WcPolicy;
#else
        using TimingPolicy_T = walberla::timing::DeviceSynchronizePolicy;
#endif

        walberla::timeloop::SweepTimeloop<TimingPolicy_T> precursorTimeLoop( precursorBlocks->getBlockStorage(), timesteps );
        walberla::timeloop::SweepTimeloop<TimingPolicy_T> concurrentTimeLoop( concurrentBlocks->getBlockStorage(), timesteps );
        farm.setTimeloop(&concurrentTimeLoop);

        // driving forces
        walberla::wind::FlowDriverCollection precursorFlowDriver(precursorBlocks, &precursorTimeLoop, globalConfig, precursorDomainSetup, precursorForceFieldGpuID, precursorVelocityFieldGpuID, fieldGhostLayers);
        for( auto & block : *precursorBlocks) {
            precursorFlowDriver(&block);
        }
        walberla::wind::FlowDriverCollection concurrentFlowDriver(concurrentBlocks, &concurrentTimeLoop, globalConfig, concurrentDomainSetup, concurrentForceFieldGpuID, concurrentVelocityFieldGpuID, fieldGhostLayers);
        for( auto & block : *concurrentBlocks) {
            concurrentFlowDriver(&block);
        }


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
        auto syncForces = [&](){farm.syncNextNeighbour(mpi::SYNC_FORCE, maxKernelWidth+uint_t(real_t(diameter)/real_t(2.0)), cudaEnabledMPI);};
        auto syncMacros = [&](){farm.syncNextNeighbour(mpi::SYNC_MACRO, 1, cudaEnabledMPI);};
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
        //eddy viscosity and strainRate computation will be done on the CPU side
        walberla::pystencils::waLBerlaWind_StrainRateWriter strainRateWriter(concurrentForceFieldGpuID, concurrentNutFieldGpuID, concurrentPdfFieldGpuID, concurrentStrainRateFieldGpuID, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2]);
        walberla::pystencils::waLBerlaWind_EddyViscosityWriter eddyViscosityWriter(concurrentNutFieldGpuID, gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2], omega);

        // getter sweep for output
        walberla::pystencils::waLBerlaWind_MacroGetter precursorGetter(precursorDensityFieldCpuID, precursorForceFieldCpuID, precursorPdfFieldCpuID, precursorVelocityFieldCpuID);
        walberla::pystencils::waLBerlaWind_MacroGetter concurrentGetter(concurrentDensityFieldCpuID, concurrentForceFieldCpuID, concurrentPdfFieldCpuID, concurrentVelocityFieldCpuID);

        //Copy the data from GPU to CPU
        concurrentTurbineOutput.addMacroscopicFieldCalculator([&]() {
            walberla::gpu::fieldCpy<ScalarField_T , GPUField_T<real_t> >( concurrentBlocks, concurrentDensityFieldCpuID, concurrentDensityFieldGpuID );
            walberla::gpu::fieldCpy<VectorField_T , GPUField_T<real_t> >( concurrentBlocks, concurrentVelocityFieldCpuID, concurrentVelocityFieldGpuID );
            walberla::gpu::fieldCpy<VectorField_T , GPUField_T<real_t> >( concurrentBlocks, concurrentForceFieldCpuID, concurrentForceFieldGpuID );

            walberla::gpu::fieldCpy<VectorField_T , GPUField_T<real_t> >( concurrentBlocks, concurrentMeanVelocityFieldWFBCpuID, concurrentMeanVelocityFieldWFBGpuID );
            walberla::gpu::fieldCpy<VectorField_T , GPUField_T<real_t> >( concurrentBlocks, concurrentMeanVelocityFieldOutputCpuID, concurrentMeanVelocityFieldOutputGpuID );
            walberla::gpu::fieldCpy<SecondOrderTensorField_T , GPUField_T<real_t> >(concurrentBlocks, concurrentSumOfSquaresFieldCpuID, concurrentSumOfSquaresFieldGpuID );

            cudaDeviceSynchronize();
        });

        precursorTurbineOutput.addMacroscopicFieldCalculator([&]() {
            walberla::gpu::fieldCpy<ScalarField_T , GPUField_T<real_t> >( precursorBlocks, precursorDensityFieldCpuID, precursorDensityFieldGpuID );
            walberla::gpu::fieldCpy<VectorField_T , GPUField_T<real_t> >( precursorBlocks, precursorVelocityFieldCpuID, precursorVelocityFieldGpuID );
            walberla::gpu::fieldCpy<VectorField_T , GPUField_T<real_t> >( precursorBlocks, precursorForceFieldCpuID, precursorForceFieldGpuID );

            walberla::gpu::fieldCpy<VectorField_T , GPUField_T<real_t> >( precursorBlocks, precursorMeanVelocityFieldWFBCpuID, precursorMeanVelocityFieldWFBGpuID );

            cudaDeviceSynchronize();
        });

        auto concurrentOutput = [&concurrentTurbineOutput](){cudaDeviceSynchronize(); concurrentTurbineOutput.write();};
        auto precursorOutput = [&precursorTurbineOutput](){cudaDeviceSynchronize(); precursorTurbineOutput.write();};

        //NOTE must convert sweeps that are altered to lambdas, otherwise copy and counter will stay 0
        precursorWelfordSweepWFB.setCounter(real_t(initialWelfordCounter));
        concurrentWelfordSweepWFB.setCounter(real_t(initialWelfordCounter));
        concurrentWelfordSweepOutput.setCounter(real_t(initialWelfordCounter));
        auto welfordLambdaWFB = [&precursorWelfordSweepWFB](walberla::IBlock * block) {
            precursorWelfordSweepWFB(block);
        };
        auto concurrentWelfordLambdaWFB = [&concurrentWelfordSweepWFB](walberla::IBlock * block) {
            concurrentWelfordSweepWFB(block);
        };
        auto welfordLambdaOutput = [&concurrentWelfordSweepOutput](walberla::IBlock * block) {
            concurrentWelfordSweepOutput(block);
        };

        concurrentTimeLoop.add() << walberla::Sweep( concurrentBoundaryCollection.getSweep(), "Concurrent boundary handling" )
        << walberla::AfterFunction( [&]() {
            auto precursorBlock = precursorBlocks->begin();
             for (auto concurrentBlock = concurrentBlocks->begin(); concurrentBlock != concurrentBlocks->end();
             ++concurrentBlock, ++precursorBlock) {

                 auto concurrentPDFField = concurrentBlock->getData< walberla::gpu::GPUField<real_t> >(concurrentPdfFieldGpuID);
                 auto precursorPDFField  = precursorBlock->getData< walberla::gpu::GPUField<real_t> >(precursorPdfFieldGpuID);
                 cudaDeviceSynchronize();

                int threadsPerBlock = 256;
                 auto nInflowCells = inflowCellsIndices[&concurrentBlock->getId()].size();
                 int blocksPerGrid = (nInflowCells+threadsPerBlock-1) / threadsPerBlock;

                copyBoundaryData<<< blocksPerGrid, threadsPerBlock>>>( (real_t*)precursorPDFField->data(), (real_t*)concurrentPDFField->data(),
                    gpuInflowCellsIndices[&concurrentBlock->getId()], nInflowCells, Stencil_T::Q, concurrentPDFField->fStride() );
             }
        });
        concurrentTimeLoop.add() << walberla::Sweep( concurrentSweepCollection.streamCollide() , "Concurrent LBM streamCollide sweep")
                       << walberla::AfterFunction( concurrentCommunication.getCommunicateFunctor(), "field communication" )
                        << walberla::AfterFunction( concurrentOutput, "Concurrent VTK output" );

        concurrentTimeLoop.add()
            << walberla::BeforeFunction([&](){
                              if(welfordIntervalOutput // && int(welfordSweepOutput.counter_)
                                 && (uint_t(concurrentWelfordSweepOutput.getCounter()) % welfordIntervalOutput == 0)){
                                  concurrentWelfordSweepOutput.setCounter(real_t(1));
                                  // reset the field: swap pointers

                                  // Take care this is called only after the VTK output, otherwise we write the wrong data in the meanVelocity output
                                  // We live a dangerous life: we swap the velocity field counter
                                  for( auto blockIt = concurrentBlocks->begin(); blockIt != concurrentBlocks->end(); ++blockIt ) {
                                      auto velocityField = blockIt->getData<walberla::gpu::GPUField<real_t>>( concurrentVelocityFieldGpuID );
                                      auto meanVelocityField = blockIt->getData<walberla::gpu::GPUField<real_t>>( concurrentMeanVelocityFieldOutputGpuID );
                                      velocityField->swapDataPointers(meanVelocityField);

                                      sopResetter(blockIt.get());
                                  }
                              } else {
                                  concurrentWelfordSweepOutput.setCounter(real_t(concurrentWelfordSweepOutput.getCounter()+1));
                              }
                          }, "welford sweep Output")
                       << walberla::Sweep(welfordLambdaOutput, "welford sweep Output");

        concurrentTimeLoop.add()
                        << walberla::BeforeFunction( [&]() {
                              farm.callback(Component::Function::UPDATE_DISCRETISATION);
                          }, "Rotation" )
                       << walberla::Sweep( evaluateDensityAndVelocity, "Evaluate density and velocity", sets::turbineSelector() )
                       << walberla::AfterFunction(syncMacros, "Communicate macroscopic turbine data")
                       << walberla::AfterFunction( applyControl , "Apply Control");

        concurrentTimeLoop.add() << walberla::Sweep( concurrentFlowDriver, "Concurrent setting driving force" );
        concurrentTimeLoop.add() << walberla::Sweep(concurrentWelfordLambdaWFB, "Concurrent welford sweep WFB")
                                 << walberla::BeforeFunction([&](){
                                       concurrentWelfordSweepWFB.setCounter(real_t(concurrentWelfordSweepWFB.getCounter()+1));
                                   }, "Concurrent welford sweep WFB");
        concurrentTimeLoop.add()
                       << walberla::BeforeFunction( [&]() {
                              farm.callback( Component::Function::CALCULATE_FORCES );
                          }, "Calculate forces" )
                       << walberla::BeforeFunction( syncForces, "Communicate force data" )
                       << walberla::Sweep(spreadForces, "Spread forces", sets::turbineSelector() )
                       << walberla::AfterFunction( [&]() {
                           farm.callback(Component::Function::RESET_COMMUNICATION_DATA);
                       }, "ResetCommunicationData" );

        // LBM stability check
        concurrentTimeLoop.addFuncAfterTimeStep( walberla::makeSharedFunctor( walberla::field::makeStabilityChecker< PdfField_T, FlagField_T >(globalConfig, concurrentBlocks, concurrentPdfFieldCpuID, concurrentFlagFieldCpuID, FluidFlagUID ) ),
                                       "Concurrent LBM stability check" );




        //precursorTimeLoop.add()
        //                << walberla::Sweep( precursorBoundaryCollection.getSweep(), "Precursor boundary handling" );
        precursorTimeLoop.add() << walberla::BeforeFunction( precursorCommunication.getCommunicateFunctor(), "field communication" )
                       << walberla::BeforeFunction( [&precursorBoundarySetup, &shiftedPeriodicity]() {
                                           if( precursorBoundarySetup.inflowType() == InflowSetup::ShiftedPeriodic )
                                               shiftedPeriodicity();
                                      }, "Shifted periodicity" )
                       << walberla::Sweep( precursorBoundaryCollection.getSweep(), "Boundary handling" );
        precursorTimeLoop.add() << walberla::Sweep( precursorSweepCollection.streamCollide() , "Precursor LBM streamCollide sweep")
                        << walberla::AfterFunction( precursorOutput, "Precursor VTK output" )

                        << walberla::AfterFunction( [&](){

                            if (precursorTimeLoop.getCurrentTimeStep() % precursorDomainSetup.snapshotHandler_->snapshotFrequency == 0 && precursorTimeLoop.getCurrentTimeStep() > 0) {
                                walberla::gpu::fieldCpy < PdfField_T, walberla::gpu::GPUField<real_t> >
                                        (precursorBlocks, precursorPdfFieldCpuID, precursorPdfFieldGpuID);
                                walberla::gpu::fieldCpy < VectorField_T, walberla::gpu::GPUField<real_t> >
                                        (precursorBlocks, precursorMeanVelocityFieldWFBCpuID, precursorMeanVelocityFieldWFBGpuID);
                                walberla::gpu::fieldCpy < VectorField_T, walberla::gpu::GPUField<real_t> >
                                        (precursorBlocks, precursorForceFieldCpuID, precursorForceFieldGpuID);
                                cudaDeviceSynchronize();
                                precursorDomainSetup.snapshotHandler_->createSnapshot(precursorTimeLoop.getCurrentTimeStep(),
                                    precursorBlocks, precursorPdfFieldCpuID, precursorForceFieldCpuID,
                        uint_t(precursorWelfordSweepWFB.getCounter()), precursorMeanVelocityFieldWFBCpuID);
                            }
                        }, "Snapshot creation");

        precursorTimeLoop.add() << walberla::BeforeFunction([&](){
                              precursorWelfordSweepWFB.setCounter(real_t(precursorWelfordSweepWFB.getCounter()+1));
                          }, "Precursor welford sweep WFB")
                       << walberla::Sweep(welfordLambdaWFB, "Precursor welford sweep WFB");
        precursorTimeLoop.add() << walberla::Sweep( precursorFlowDriver, "Precursor setting driving force" );

        precursorTimeLoop.addFuncAfterTimeStep( walberla::timing::RemainingTimeLogger( precursorTimeLoop.getNrOfTimeSteps(), remainingTimeLoggerFrequency ), "Precursor remaining time logger" );

        precursorTimeLoop.addFuncAfterTimeStep( walberla::makeSharedFunctor( walberla::field::makeStabilityChecker< PdfField_T, FlagField_T >(globalConfig, precursorBlocks, precursorPdfFieldCpuID, precursorFlagFieldCpuID, FluidFlagUID ) ),
                                       "Precursor LBM stability check" );

        WALBERLA_LOG_INFO_ON_ROOT("Running timeloop...")

        WALBERLA_MPI_WORLD_BARRIER()

        walberla::timing::TimingPool<TimingPolicy_T> concurrentTiming;
        walberla::timing::TimingPool<TimingPolicy_T> precursorTiming;
        concurrentTiming.registerTimer("Force Output");

        walberla::WcTimer concurrentTimer;
        walberla::WcTimer precursorTimer;
        precursorTimer.start();

        for( uint_t i = 0; i < timesteps; ++i ) {

            precursorTimeLoop.singleStep(precursorTiming);

            if(i == concurrentStart) {
                WALBERLA_LOG_INFO_ON_ROOT("Initialising concurrent simulation...")

                concurrentTimer.start();

                walberla::gpu::fieldCpy<PdfField_T, walberla::gpu::GPUField<real_t> >( precursorBlocks, precursorPdfFieldCpuID, precursorPdfFieldGpuID );
                walberla::gpu::fieldCpy<ScalarField_T, walberla::gpu::GPUField<real_t> >( precursorBlocks, precursorDensityFieldCpuID, precursorDensityFieldGpuID );
                walberla::gpu::fieldCpy<VectorField_T, walberla::gpu::GPUField<real_t> >( precursorBlocks, precursorVelocityFieldCpuID, precursorVelocityFieldGpuID );
                walberla::gpu::fieldCpy<VectorField_T, walberla::gpu::GPUField<real_t> >( precursorBlocks, precursorMeanVelocityFieldWFBCpuID, precursorMeanVelocityFieldWFBGpuID );
                walberla::gpu::fieldCpy<VectorField_T, walberla::gpu::GPUField<real_t> >( precursorBlocks, precursorForceFieldCpuID, precursorForceFieldGpuID );
                auto concurrentBlock = concurrentBlocks->begin();
                for (auto precursorBlock = precursorBlocks->begin(); precursorBlock != precursorBlocks->end();
                        ++precursorBlock) {

                        auto precursorPDFField  = precursorBlock->getData< PdfField_T >(precursorPdfFieldCpuID);
                        auto concurrentPDFField   = concurrentBlock->getData< PdfField_T >(concurrentPdfFieldCpuID);

                        auto precursorDensityField  = precursorBlock->getData< ScalarField_T >(precursorDensityFieldCpuID);
                        auto concurrentDensityField   = concurrentBlock->getData< ScalarField_T >(concurrentDensityFieldCpuID);

                        auto precursorVelocityField  = precursorBlock->getData< VectorField_T >(precursorVelocityFieldCpuID);
                        auto concurrentVelocityField   = concurrentBlock->getData< VectorField_T >(concurrentVelocityFieldCpuID);

                        auto precursorWFBField = precursorBlock->getData< VectorField_T >(precursorMeanVelocityFieldWFBCpuID);
                        auto concurrentWFBField   = concurrentBlock->getData< VectorField_T >(concurrentMeanVelocityFieldWFBCpuID);

                        auto precursorForceField  = precursorBlock->getData< VectorField_T >(precursorForceFieldCpuID);
                        auto concurrentForceField   = concurrentBlock->getData< VectorField_T >(concurrentForceFieldCpuID);

                        auto concurrentCellIt = concurrentPDFField->beginWithGhostLayer();
                        for (auto precursorCellIt = precursorPDFField->beginWithGhostLayer();
                                precursorCellIt != precursorPDFField->end(); ++precursorCellIt)
                        {
                            const walberla::Cell localPrecursorCell = precursorCellIt.cell();
                            const walberla::Cell localConcurrentCell  = concurrentCellIt.cell();

                            for (uint_t idx = 0; idx < Stencil_T::Q; idx++)
                            {
                                concurrentPDFField->get(localConcurrentCell, idx) =
                                        precursorPDFField->get(localPrecursorCell, idx);
                            }
                            concurrentDensityField->get(localConcurrentCell) =
                                precursorDensityField->get(localPrecursorCell);
                            for (uint_t idx = 0; idx < 3; idx++)
                            {
                                concurrentVelocityField->get(localConcurrentCell, idx) =
                                        precursorVelocityField->get(localPrecursorCell, idx);
                                concurrentWFBField->get(localConcurrentCell, idx) =
                                        precursorWFBField->get(localPrecursorCell, idx);
                                concurrentForceField->get(localConcurrentCell, idx) =
                                        precursorForceField->get(localPrecursorCell, idx);
                            }
                            ++concurrentCellIt;
                        }
                        concurrentBlock++;
                }
                walberla::gpu::fieldCpy<walberla::gpu::GPUField<real_t>, PdfField_T>( concurrentBlocks, concurrentPdfFieldGpuID, concurrentPdfFieldCpuID );
                walberla::gpu::fieldCpy<walberla::gpu::GPUField<real_t>, ScalarField_T>( concurrentBlocks, concurrentDensityFieldGpuID, concurrentDensityFieldCpuID );
                walberla::gpu::fieldCpy<walberla::gpu::GPUField<real_t>, VectorField_T>( concurrentBlocks, concurrentVelocityFieldGpuID, concurrentVelocityFieldCpuID );
                walberla::gpu::fieldCpy<walberla::gpu::GPUField<real_t>, VectorField_T>( concurrentBlocks, concurrentMeanVelocityFieldWFBGpuID, concurrentMeanVelocityFieldWFBCpuID );
                walberla::gpu::fieldCpy<walberla::gpu::GPUField<real_t>, VectorField_T>( concurrentBlocks, concurrentForceFieldGpuID, concurrentForceFieldCpuID );
                cudaDeviceSynchronize();
            }

            if (i >= concurrentStart)
            {
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

        walberla::lbm::PerformanceEvaluation<FlagField_T> concurrentPerformance(concurrentBlocks, concurrentFlagFieldCpuID, FluidFlagUID);
        concurrentPerformance.logResultOnRoot(timesteps, concurrentTime);
        walberla::lbm::PerformanceEvaluation<FlagField_T> precursorPerformance(precursorBlocks, precursorFlagFieldCpuID, FluidFlagUID);
        precursorPerformance.logResultOnRoot(timesteps, precursorTime);

        farm.calculatePowerAndThrust();

        logging::WindTurbineApplicationLogger<domain::DomainSetup, domain::BoundarySetup, codegen::KernelInfo, TimingPolicy_T> concurrentLogger {
            concurrentDomainSetup, concurrentBoundarySetup, windFarmConfig, &concurrentTiming, concurrentPerformance.loggingString(timesteps, concurrentTime),
            cudaEnabledMPI, "WindTurbine_Logging_concurrent.log"
        };
        logging::WindTurbineApplicationLogger<domain::DomainSetup, domain::BoundarySetup, codegen::KernelInfo, TimingPolicy_T> precursorLogger {
            precursorDomainSetup, precursorBoundarySetup, windFarmConfig, &precursorTiming, precursorPerformance.loggingString(timesteps, precursorTime),
            cudaEnabledMPI, "WindTurbine_Logging_precursor.log"
        };

         for (auto concurrentBlock = concurrentBlocks->begin(); concurrentBlock != concurrentBlocks->end();
              ++concurrentBlock) {
             cudaFree(gpuInflowCellsIndices[&(concurrentBlock->getId())]);
        }

        return EXIT_SUCCESS;

    }

}

int main(int argc, char** argv) {
    return turbine_core::main(argc, argv);
}
