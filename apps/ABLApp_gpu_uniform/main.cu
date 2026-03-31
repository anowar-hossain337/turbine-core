#include <iostream>
#include <cmath>
#include <vector>

#include <blockforest/Initialization.h>
#include <core/all.h>
#include <core/Environment.h>
#include <core/mpi/all.h>
#include <domain_decomposition/all.h>
#include <field/all.h>
#include <geometry/all.h>
#include <lbm/all.h>
#include <lbm_generated/field/AddToStorage.h>
#include <timeloop/all.h>
#include <vtk/VTKOutput.h>  // VTK output (GPU-safe): copy selected fields to CPU before writing
#include <field/vtk/VTKWriter.h>  // VTK output (GPU-safe): copy selected fields to CPU before writing

#include "conversion/Conversion.h"
#include "domain/BoundaryHandling.h"
#include "domain/BoundarySetup.h"
#include "domain/DomainInitialiser.h"
#include "domain/DomainSetup.h"
#include "walberla_helper/field/all.h"
#include "wind_turbine_core/ProjectDefines.h"

#include "waLBerlaABL_KernelInfo.h"
#include "FlowDriverCollection.h"

namespace turbine_core {

static const uint_t fieldGhostLayers = 1;

template<typename Type_T>
using GPUField_T = walberla::gpu::GPUField<Type_T>;

int main(int argc, char** argv) {

    walberla::Environment walberlaEnv(argc, argv);
    walberla::logging::Logging::instance()->setLogLevel(walberla::logging::Logging::INFO);

    auto globalConfig = walberlaEnv.config();
    auto parameters = globalConfig->getOneBlock("Parameters");

    uint_t timesteps = parameters.getParameter<uint_t>("timesteps", uint_t(10));
    ++timesteps;

    const real_t remainingTimeLoggerFrequency = parameters.getParameter<real_t>("remainingTimeLoggerFrequency", real_t(3.0));
    const real_t length_SI = parameters.getParameter<real_t>("length_SI");
    const real_t length_LU = parameters.getParameter<real_t>("length_LU");
    const real_t velocity_SI = parameters.getParameter<real_t>("velocity_SI");
    const real_t velocity_LU = parameters.getParameter<real_t>("velocity_LU");
    const real_t viscosity_SI = parameters.getParameter<real_t>("viscosity_SI");
    const real_t density_SI = parameters.getParameter<real_t>("density_SI");

    Conversion::calculateConversionFactors(length_SI, length_LU, velocity_SI, velocity_LU, density_SI);
    Conversion::print();

    const real_t viscosity_LU = (viscosity_SI / density_SI) * Conversion::C_t() / Conversion::C_l() / Conversion::C_l();
    WALBERLA_CHECK_GREATER(viscosity_LU, real_t(0.0), "Negative LU velocity - Did you calculate the conversion factors?")

    const real_t omega = walberla::lbm::collision_model::omegaFromViscosity(viscosity_LU);
    WALBERLA_LOG_INFO_ON_ROOT("omega = " << omega)

    WALBERLA_LOG_INFO_ON_ROOT("Creating blockforest...")

    auto boundariesConfig = globalConfig->getOneBlock("Boundaries");
    domain::BoundarySetup boundarySetup(boundariesConfig);

    domain::DomainSetup domainSetup(globalConfig, boundarySetup.periodicity());

    WALBERLA_CHECK(!domainSetup.snapshotHandler_->loadSnapshot,
                   "ABLApp_gpu_uniform does not support snapshot restart. Please set DomainSetup::loadSnapshot = 0.")
    WALBERLA_CHECK(!domainSetup.snapshotHandler_->storeSnapshot,
                   "ABLApp_gpu_uniform does not support snapshot writing. Please set DomainSetup::storeSnapshot = 0.")

    std::vector<AABB> emptyAABBs;
    auto blocks = domainSetup.createUniformBlockForest(emptyAABBs);

    WALBERLA_LOG_INFO_ON_ROOT("Creating fields...")

    auto layout = codegen::KernelInfo::layout;
    const StorageSpecification_T storageSpec{};

    auto allocator = std::make_shared<walberla::gpu::HostFieldAllocator<real_t>>();

    const BlockDataID forceFieldCpuID = walberla::field::addToStorage<VectorField_T>(blocks, "force field CPU", real_t(0), layout, fieldGhostLayers, allocator);
    const BlockDataID meanVelocityWfbFieldCpuID = walberla::field::addToStorage<VectorField_T>(blocks, "mean velocity wfb field CPU", real_t(0), layout, fieldGhostLayers, allocator);
    const BlockDataID pdfFieldCpuID = walberla::lbm_generated::addPdfFieldToStorage(
            blocks, "pdf field CPU", storageSpec, fieldGhostLayers, layout,
            walberla::Set<walberla::SUID>::emptySet(), walberla::Set<walberla::SUID>::emptySet(), allocator);

    const BlockDataID densityFieldCpuID = walberla::field::addToStorage<ScalarField_T>(blocks, "density field CPU", real_t(1), layout, fieldGhostLayers, allocator);
    const BlockDataID eddyViscosityFieldCpuID = walberla::field::addToStorage<ScalarField_T>(blocks, "eddy viscosity field CPU", real_t(0), layout, fieldGhostLayers, allocator);
    const BlockDataID flagFieldID = walberla::field::addFlagFieldToStorage<FlagField_T>(blocks, "flag field", fieldGhostLayers);
    const BlockDataID omegaFieldCpuID = walberla::field::addToStorage<ScalarField_T>(blocks, "omega field CPU", real_t(0), layout, fieldGhostLayers, allocator);
    const BlockDataID velocityFieldCpuID = walberla::field::addToStorage<VectorField_T>(blocks, "velocity field CPU", real_t(0), layout, fieldGhostLayers, allocator);

    auto initParams = globalConfig->getOneBlock("Initialisation");
    const auto initType = domain::DomainInitialisation::toType(initParams.getParameter<std::string>("type"));
    const Vector3<real_t> initialVelocity = initParams.getParameter<Vector3<real_t>>("initialVelocity");

    const real_t kappa = parameters.getParameter<real_t>("kappa", real_t(0.42));
    const real_t B = parameters.getParameter<real_t>("B", real_t(5.5));
    const real_t roughnessLengthRatio = parameters.getParameter<real_t>("roughnessLengthRatio", real_t(1e-4));
    const real_t referenceHeight = parameters.getParameter<real_t>("referenceHeight_LU", real_t(-1));
    const uint32_t samplingHeight = parameters.getParameter<uint32_t>("samplingHeight_LU", uint32_t(0));
    const real_t roughnessLength = roughnessLengthRatio * referenceHeight;

    const real_t uTau = kappa * initialVelocity[0] / std::log(real_t(1) / roughnessLengthRatio);

    WALBERLA_LOG_INFO_ON_ROOT("Initialise domain...")

    std::unique_ptr<domain::DomainInitialisation> initialiser{};
    if(initType == domain::DomainInitialisation::UNIFORM) {
        initialiser = std::make_unique<domain::UniformInitialisation>(blocks, initialVelocity);
    } else if(initType == domain::DomainInitialisation::ASMUTH) {
        initialiser = std::make_unique<domain::AsmuthInitialisation>(blocks, uTau, roughnessLength, kappa, B,
                                                                     viscosity_LU, domainSetup.domainSize_, true);
    } else if(initType == domain::DomainInitialisation::LOG_LAW) {
        initialiser = std::make_unique<domain::AsmuthInitialisation>(blocks, uTau, roughnessLength, kappa, B,
                                                                     viscosity_LU, domainSetup.domainSize_, false);
    } else {
        WALBERLA_ABORT("Unsupported initialisation type")
    }

    using PdfSetter_T = walberla::pystencils::waLBerlaABL_PdfSetter;
    PdfSetter_T setter(densityFieldCpuID, forceFieldCpuID, pdfFieldCpuID, velocityFieldCpuID);
    initialiser->setViaVelocityField<VectorField_T, PdfSetter_T>(velocityFieldCpuID, setter);
    initialiser->setViaVelocityField<VectorField_T, PdfSetter_T>(meanVelocityWfbFieldCpuID, setter);

    WALBERLA_LOG_INFO_ON_ROOT("Add GPU fields...")

    const BlockDataID densityFieldGpuID = walberla::gpu::addGPUFieldToStorage<ScalarField_T>(blocks, densityFieldCpuID, "density field GPU", true);
    const BlockDataID eddyViscosityFieldGpuID = walberla::gpu::addGPUFieldToStorage<ScalarField_T>(blocks, eddyViscosityFieldCpuID, "eddy viscosity field GPU", true);
    const BlockDataID forceFieldGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T>(blocks, forceFieldCpuID, "force field GPU", true);
    const BlockDataID meanVelocityWfbFieldGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T>(blocks, meanVelocityWfbFieldCpuID, "mean velocity wfb field GPU", true);
    const BlockDataID omegaFieldGpuID = walberla::gpu::addGPUFieldToStorage<ScalarField_T>(blocks, omegaFieldCpuID, "omega field GPU", true);
    const BlockDataID pdfFieldGpuID = walberla::lbm_generated::addGPUPdfFieldToStorage<PdfField_T, StorageSpecification_T>(blocks, pdfFieldCpuID, storageSpec, "pdf field GPU");
    const BlockDataID velocityFieldGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T>(blocks, velocityFieldCpuID, "velocity field GPU", true);

    SweepCollection_T sweepCollection(blocks, densityFieldGpuID, eddyViscosityFieldGpuID, forceFieldGpuID, omegaFieldGpuID,
                                      pdfFieldGpuID, velocityFieldGpuID, omega);

    WALBERLA_MPI_BARRIER()
    WALBERLA_LOG_INFO_ON_ROOT("Initialisation done")

    WALBERLA_LOG_INFO_ON_ROOT("Setting up boundaries...")

    boundarySetup.fillFlagFieldFromConfig<FlagField_T>(blocks, flagFieldID, FluidFlagUID,
                                                        NoSlipFlagUID, WFBFlagUID, SymmetryFlagUID, UniformInflowFlagUID,
                                                        LogLawInflowFlagUID, OutflowFlagUID);

    if(boundarySetup.wallType() == WallSetup::WFB) {
        WALBERLA_CHECK(referenceHeight > real_t(0), "Reference height must be given in the parameter file for wall function boundary conditions.")
        WALBERLA_CHECK(boundarySetup.environmentSetup() != EnvironmentSetup::Tunnel,
                       "Wall-function bounce is currently not supported for the tunnel environment")
        WALBERLA_CHECK(samplingHeight > uint_t(0), "Sampling height must be given in the parameter file for wall function boundary conditions.")
    }

    const auto& inflowVelocity = boundarySetup.inflowVelocity();
    auto velocityInit = boundary::velocityInit(roughnessLength, kappa, uTau);
    BoundaryCollection_T boundaryCollection(
            blocks, flagFieldID, pdfFieldGpuID, FluidFlagUID, forceFieldGpuID,
            meanVelocityWfbFieldGpuID, samplingHeight,
            roughnessLengthRatio * referenceHeight,
            inflowVelocity[0], inflowVelocity[1], inflowVelocity[2], velocityInit,
            pdfFieldCpuID
    );

    const int periodicShiftValue = boundariesConfig.getParameter<int>(
            "periodicShiftValue", int(real_t(0.33) * real_t(domainSetup.domainSize_[2])));
    walberla::gpu::ShiftedPeriodicityGPU<PdfGPUField_T> shiftedPeriodicity(blocks, pdfFieldGpuID, fieldGhostLayers,
                                                                           0, 1, periodicShiftValue);

    walberla::pystencils::waLBerlaABL_WelfordWFB welfordWFBSweep(
            meanVelocityWfbFieldGpuID, velocityFieldGpuID,
            real_t(0));
    auto welfordWFBLambda = [&welfordWFBSweep](walberla::IBlock * block) { welfordWFBSweep(block); };

    WALBERLA_LOG_INFO_ON_ROOT("Set up communication...")

    bool cudaEnabledMPI = parameters.getParameter<bool>("cudaEnabledMPI", false);
    if(walberla::mpi::MPIManager::instance()->numProcesses() == 1) {
        cudaEnabledMPI = false;
    }

    auto communication = std::make_shared<walberla::gpu::communication::UniformGPUScheme<CommunicationStencil_T>>(blocks, cudaEnabledMPI);
    auto pdfPackInfo = std::make_shared<walberla::lbm_generated::UniformGeneratedGPUPdfPackInfo<PdfGPUField_T>>(pdfFieldGpuID);
    communication->addPackInfo(pdfPackInfo);
    auto densityPackInfo = std::make_shared<walberla::gpu::communication::MemcpyPackInfo<GPUField_T<real_t>>>(densityFieldGpuID);
    communication->addPackInfo(densityPackInfo);
    auto velocityPackInfo = std::make_shared<walberla::gpu::communication::MemcpyPackInfo<GPUField_T<real_t>>>(velocityFieldGpuID);
    communication->addPackInfo(velocityPackInfo);

    WALBERLA_LOG_INFO_ON_ROOT("Creating time loop...")

#ifdef NDEBUG
    using TimingPolicy_T = walberla::timing::WcPolicy;
#else
    using TimingPolicy_T = walberla::timing::DeviceSynchronizePolicy;
#endif

    auto timeloop = walberla::timeloop::SweepTimeloop<TimingPolicy_T>(blocks->getBlockStorage(), timesteps);

    // VTK output (GPU-safe): copy selected fields to CPU before writing
    auto outputConfig = globalConfig->getOneBlock("Output");
    auto vtkConfig = outputConfig.getOneBlock("VTK");
    const uint_t vtkWriteFrequency = vtkConfig.getParameter<uint_t>("writeFrequency", uint_t(0));
    const uint_t vtkStartTimestep = vtkConfig.getParameter<uint_t>("startTimestep", uint_t(0));
    const uint_t vtkGhostLayers = vtkConfig.getParameter<uint_t>("ghostLayers", uint_t(0));
    const std::string vtkBaseFolder = vtkConfig.getParameter<std::string>("baseFolder", "vtk_out");
    const std::string vtkExecutionFolder = vtkConfig.getParameter<std::string>("executionFolder", "simulation_step");

    std::shared_ptr<walberla::vtk::VTKOutput> fieldVTKOutput{nullptr};
    if(vtkWriteFrequency > 0) {
        fieldVTKOutput = walberla::vtk::createVTKOutput_BlockData(
                *blocks,
                "abl_gpu_uniform",
                vtkWriteFrequency,
                vtkGhostLayers,
                false,
                vtkBaseFolder,
                vtkExecutionFolder,
                true,
                true,
                true,
                true,
                uint_t(0),
                false,
                false);

        // fieldVTKOutput->addCellDataWriter(std::make_shared<walberla::field::VTKWriter<ScalarField_T>>(densityFieldCpuID, "Density"));
        fieldVTKOutput->addCellDataWriter(std::make_shared<walberla::field::VTKWriter<VectorField_T>>(velocityFieldCpuID, "Velocity"));
        // fieldVTKOutput->addCellDataWriter(std::make_shared<walberla::field::VTKWriter<VectorField_T>>(forceFieldCpuID, "Force"));
        // fieldVTKOutput->addCellDataWriter(std::make_shared<walberla::field::VTKWriter<ScalarField_T>>(eddyViscosityFieldCpuID, "EddyViscosity"));
        // fieldVTKOutput->addCellDataWriter(std::make_shared<walberla::field::VTKWriter<ScalarField_T>>(omegaFieldCpuID, "Omega"));
        // fieldVTKOutput->addCellDataWriter(std::make_shared<walberla::field::VTKWriter<FlagField_T>>(flagFieldID, "Flag"));
    }

    auto writeVTK = [&]() {
        if(!fieldVTKOutput) {
            return;
        }

        static uint_t vtkStepCounter = uint_t(0);
        if(vtkStepCounter < vtkStartTimestep) {
            ++vtkStepCounter;
            return;
        }

        //walberla::gpu::fieldCpy<ScalarField_T, GPUField_T<real_t>>(blocks, densityFieldCpuID, densityFieldGpuID);
        walberla::gpu::fieldCpy<VectorField_T, GPUField_T<real_t>>(blocks, velocityFieldCpuID, velocityFieldGpuID);
        //walberla::gpu::fieldCpy<VectorField_T, GPUField_T<real_t>>(blocks, forceFieldCpuID, forceFieldGpuID);
        //walberla::gpu::fieldCpy<ScalarField_T, GPUField_T<real_t>>(blocks, eddyViscosityFieldCpuID, eddyViscosityFieldGpuID);
        //walberla::gpu::fieldCpy<ScalarField_T, GPUField_T<real_t>>(blocks, omegaFieldCpuID, omegaFieldGpuID);
        cudaDeviceSynchronize();

        fieldVTKOutput->write();
        ++vtkStepCounter;
    };
// VTK output (GPU-safe): copy selected fields to CPU before writing
    walberla::wind::FlowDriverCollection flowDriver(blocks, &timeloop, globalConfig, domainSetup, forceFieldGpuID, velocityFieldGpuID, fieldGhostLayers);
    for(auto & block : *blocks) { flowDriver(&block); }

    timeloop.add() << walberla::BeforeFunction(communication->getCommunicateFunctor(), "Field communication")
                   << walberla::BeforeFunction([&boundarySetup, &shiftedPeriodicity]() {
                       if(boundarySetup.inflowType() == InflowSetup::ShiftedPeriodic) shiftedPeriodicity();
                   }, "Shifted periodicity")
                   << walberla::Sweep(boundaryCollection.getSweep(), "Boundary handling");

    timeloop.add() << walberla::Sweep(sweepCollection.streamCollide(), "LBM stream-collide")
                   << walberla::AfterFunction(writeVTK, "VTK output");  // VTK output (GPU-safe): copy selected fields to CPU before writing

    if(boundarySetup.wallType() == WallSetup::WFB) {
        timeloop.add() << walberla::BeforeFunction([&]() {
                           welfordWFBSweep.setCounter(real_t(welfordWFBSweep.getCounter() + 1));
                       }, "WelfordWFB counter")
                       << walberla::Sweep(welfordWFBLambda, "WelfordWFB sweep");
    }

    timeloop.add() << walberla::Sweep(flowDriver, "Setting driving force");

    timeloop.addFuncAfterTimeStep(walberla::timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
                                  "Remaining time logger");

    WALBERLA_LOG_INFO_ON_ROOT("Running timeloop...")
    WALBERLA_MPI_WORLD_BARRIER()

    walberla::timing::TimingPool<TimingPolicy_T> timing;
    walberla::WcTimer timer;
    timer.start();

    for(uint_t i = 0; i < timesteps; ++i) {
        timeloop.singleStep(timing);
    }

    timer.end();

    double time = timer.max();
    walberla::mpi::reduceInplace(time, walberla::mpi::MAX);

    const auto timeloopTiming = timing.getReduced();
    WALBERLA_LOG_INFO_ON_ROOT("Timeloop timing:\n" << *timeloopTiming)

    walberla::lbm::PerformanceEvaluation<FlagField_T> performance(blocks, flagFieldID, FluidFlagUID);
    performance.logResultOnRoot(timesteps, time);

    return EXIT_SUCCESS;
}

} // namespace turbine_core

int main(int argc, char** argv) {
    return turbine_core::main(argc, argv);
}
