// Reason for edit start: add the standard-library helpers needed to validate the new top-boundary option without changing the existing Open/WFB configuration flow.
#include <algorithm>
#include <cctype>
#include <iostream>
#include <map>
#include <cmath>
#include <string>
#include <vector>
// Reason for edit end: add the standard-library helpers needed to validate the new top-boundary option without changing the existing Open/WFB configuration flow.

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
// Reason for edit start: include the reusable filter and line-output helpers so prm-driven plane-filtered VTK and optional line sampling can be enabled without removing the existing hardcoded writer path.
#include <field/vtk/FlagFieldCellFilter.h>
#include <vtk/AABBCellFilter.h>
#include <vtk/ChainedFilter.h>
#include "wind_turbine_core/WalberlaDataTypes.h"
#include "walberla_helper/field/Field.h"
#include "output/LineOutput.h"
#include "output/PlaneInclusionFilter.h"
// Reason for edit end: include the reusable filter and line-output helpers so prm-driven plane-filtered VTK and optional line sampling can be enabled without removing the existing hardcoded writer path.
// Reason for edit start: include an app-local dynamic AABB filter so VTK inclusion boxes can move or grow without changing waLBerla's stock AABB filter.
#include "DynamicAABBInclusionFilter.h"
// Reason for edit end: include an app-local dynamic AABB filter so VTK inclusion boxes can move or grow without changing waLBerla's stock AABB filter.

#include "conversion/Conversion.h"
#include "domain/BoundaryHandling.h"
#include "domain/BoundarySetup.h"
#include "domain/DomainInitialiser.h"
#include "domain/DomainSetup.h"
#include "walberla_helper/field/all.h"
#include "wind_turbine_core/ProjectDefines.h"

// Reason for edit start: include the MesaPD and PSM runtime helpers so one optional Phase 3 kinematic sphere can be coupled on top of the existing static-obstacle ABL path when explicitly enabled.
#include "lbm_mesapd_coupling/DataTypesCodegen.h"
#include "lbm_mesapd_coupling/partially_saturated_cells_method/codegen/PSMSweepCollection.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "waLBerlaABLPSM_InitializeDomainForPSM.h"
#include "waLBerlaABLPSM_Sweep.h"
// Reason for edit end: include the MesaPD and PSM runtime helpers so one optional Phase 3 kinematic sphere can be coupled on top of the existing static-obstacle ABL path when explicitly enabled.

// Reason for edit start: include the shared moving-body config parser and both generated kernel-info headers so the app can validate and log the legacy ABL path and the optional Phase 3 PSM path from prm settings.
#include "MovingBodyConfig.h"
#include "TopDampingZone.h"
#include "TopSlipZeroGradientBoundary.h"
#include "waLBerlaABL_KernelInfo.h"
#include "waLBerlaABLPSM_KernelInfo.h"
#include "FlowDriverCollection.h"
// Reason for edit end: include the shared moving-body config parser and both generated kernel-info headers so the app can validate and log the legacy ABL path and the optional Phase 3 PSM path from prm settings.

namespace turbine_core {

static const uint_t fieldGhostLayers = 1;

template<typename Type_T>
using GPUField_T = walberla::gpu::GPUField<Type_T>;

// Reason for edit start: define the optional line-output helper once so prm-driven probe lines can be enabled while keeping the existing hardcoded VTK writer untouched.
using LineOutput_T = output::LineOutput<
        ScalarField_T, VectorField_T, PdfGPUField_T, GPUField_T<real_t>, Stencil_T,
        StorageSpecification_T::zeroCenteredPDFs, StorageSpecification_T::compressible>;
// Reason for edit end: define the optional line-output helper once so prm-driven probe lines can be enabled while keeping the existing hardcoded VTK writer untouched.

int main(int argc, char** argv) {

    walberla::Environment walberlaEnv(argc, argv);
    walberla::logging::Logging::instance()->setLogLevel(walberla::logging::Logging::INFO);

    auto globalConfig = walberlaEnv.config();
    auto parameters = globalConfig->getOneBlock("Parameters");
    // Reason for edit start: normalize prm tokens once so the optional Phase 3 moving-body checks and the existing top-boundary checks can share the same case-insensitive parsing logic.
    const auto normalizeConfigToken = [](std::string value) {
        value.erase(std::remove_if(value.begin(), value.end(),
                                   [](unsigned char ch) { return std::isspace(ch) != 0; }),
                    value.end());
        std::transform(value.begin(), value.end(), value.begin(),
                       [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
        return value;
    };
    // Reason for edit end: normalize prm tokens once so the optional Phase 3 moving-body checks and the existing top-boundary checks can share the same case-insensitive parsing logic.

    // Reason for edit start: validate the optional Phase 3 kinematic-sphere path early so the default urban ABL run stays unchanged unless all three feature toggles are enabled together.
    const moving_body::RuntimeConfig movingBodyRuntimeConfig = moving_body::loadRuntimeConfig(globalConfig);
    const bool psmCouplingRequested = movingBodyRuntimeConfig.anyFeatureEnabled();
    const bool psmCouplingEnabled =
            movingBodyRuntimeConfig.mesaPD.enabled &&
            movingBodyRuntimeConfig.psm.enabled &&
            movingBodyRuntimeConfig.body.enabled;
    const std::string movingBodyModeName = normalizeConfigToken(movingBodyRuntimeConfig.body.mode);
    const std::string movingBodyRepresentationName =
            normalizeConfigToken(movingBodyRuntimeConfig.body.representation);
    const std::string movingBodyTrajectoryTypeName =
            normalizeConfigToken(movingBodyRuntimeConfig.body.trajectoryType);
    if (psmCouplingRequested && !psmCouplingEnabled) {
        WALBERLA_ABORT("Phase 3 moving-body coupling requires MesaPD::enabled, PSM::enabled, and MovingBody::enabled to all be true together.")
    }
    if (psmCouplingEnabled) {
        WALBERLA_CHECK(!movingBodyRuntimeConfig.psm.phase1ScaffoldOnly,
                       "Phase 3 runtime coupling requires PSM::phase1ScaffoldOnly = false.")
        WALBERLA_CHECK(movingBodyModeName == "kinematic",
                       "Phase 3 currently supports only MovingBody::mode = kinematic.")
        WALBERLA_CHECK(movingBodyRepresentationName == "sphere",
                       "Phase 3 currently supports only MovingBody::representation = sphere.")
        WALBERLA_CHECK_GREATER(movingBodyRuntimeConfig.body.radius, real_t(0),
                               "MovingBody::radius must be positive.")
        if (movingBodyRuntimeConfig.body.trajectoryPointCount > uint_t(1)) {
            WALBERLA_CHECK(movingBodyTrajectoryTypeName == "piecewise_linear",
                           "Phase 3 currently supports only Trajectory::type = piecewise_linear.")
        }

        WALBERLA_LOG_INFO_ON_ROOT("MesaPD/PSM Phase 3 kinematic-sphere coupling enabled.")
        WALBERLA_LOG_INFO_ON_ROOT("PSM runtime numerics: stencil = " << codegen::psm::KernelInfo::stencil
                                  << ", method = " << codegen::psm::KernelInfo::method
                                  << ", forceModel = " << codegen::psm::KernelInfo::forceModel
                                  << ", maxParticlesPerCell = " << codegen::psm::KernelInfo::maxParticlesPerCell)

        if (movingBodyRuntimeConfig.body.trajectoryPointCount > uint_t(1)) {
            WALBERLA_LOG_INFO_ON_ROOT("Moving sphere trajectory: piecewise-linear keyframes = "
                                      << movingBodyRuntimeConfig.body.trajectoryPointCount)
        } else if (movingBodyRuntimeConfig.body.trajectoryPointCount == uint_t(1)) {
            WALBERLA_LOG_INFO_ON_ROOT("Moving sphere trajectory: single keyframe, sphere will stay fixed at "
                                      << movingBodyRuntimeConfig.body.trajectoryPoints.front().center)
        } else if (movingBodyRuntimeConfig.body.initialVelocity.sqrLength() > real_t(0)) {
            WALBERLA_LOG_INFO_ON_ROOT("Moving sphere trajectory: constant-velocity motion from initialPosition.")
        } else {
            WALBERLA_LOG_INFO_ON_ROOT("Moving sphere trajectory: no keyframes and zero initialVelocity, sphere remains stationary.")
        }
        if (movingBodyRuntimeConfig.mesaPD.subcycles != uint_t(1)) {
            WALBERLA_LOG_INFO_ON_ROOT("MesaPD::subcycles is parsed, but Phase 3 still uses one LBM update per timestep without subcycling.")
        }
    }
    // Reason for edit end: validate the optional Phase 3 kinematic-sphere path early so the default urban ABL run stays unchanged unless all three feature toggles are enabled together.

    // Reason for edit start: evaluate the optional Phase 3 kinematic sphere from either a constant velocity or piecewise-linear keyframes so the PSM path can update the particle state every timestep.
    struct MovingSphereState
    {
        walberla::Vector3<real_t> center{ real_t(0), real_t(0), real_t(0) };
        walberla::Vector3<real_t> velocity{ real_t(0), real_t(0), real_t(0) };
    };

    const auto evaluateMovingSphereState = [&](const uint_t step) {
        MovingSphereState state{};

        if (movingBodyRuntimeConfig.body.trajectoryPoints.empty()) {
            state.center = movingBodyRuntimeConfig.body.initialPosition;
            state.velocity = movingBodyRuntimeConfig.body.initialVelocity;
            for (uint_t d = uint_t(0); d < uint_t(3); ++d) {
                state.center[d] += real_t(step) * state.velocity[d];
            }
            return state;
        }

        if (movingBodyRuntimeConfig.body.trajectoryPoints.size() == size_t(1) ||
            step <= movingBodyRuntimeConfig.body.trajectoryPoints.front().step)
        {
            state.center = movingBodyRuntimeConfig.body.trajectoryPoints.front().center;
            return state;
        }

        if (step >= movingBodyRuntimeConfig.body.trajectoryPoints.back().step) {
            state.center = movingBodyRuntimeConfig.body.trajectoryPoints.back().center;
            return state;
        }

        const auto upper = std::lower_bound(
                movingBodyRuntimeConfig.body.trajectoryPoints.begin(),
                movingBodyRuntimeConfig.body.trajectoryPoints.end(),
                step,
                [](const moving_body::TrajectoryPoint & point, const uint_t testStep) {
                    return point.step < testStep;
                });

        WALBERLA_CHECK(upper != movingBodyRuntimeConfig.body.trajectoryPoints.begin(),
                       "Phase 3 trajectory interpolation requires a lower keyframe.")
        WALBERLA_CHECK(upper != movingBodyRuntimeConfig.body.trajectoryPoints.end(),
                       "Phase 3 trajectory interpolation requires an upper keyframe.")

        if (upper->step == step) {
            state.center = upper->center;
            const auto & lower = *(upper - 1);
            const real_t dt = real_t(upper->step - lower.step);
            for (uint_t d = uint_t(0); d < uint_t(3); ++d) {
                state.velocity[d] = (upper->center[d] - lower.center[d]) / dt;
            }
            return state;
        }

        const auto & lower = *(upper - 1);
        const real_t dt = real_t(upper->step - lower.step);
        const real_t alpha = real_t(step - lower.step) / dt;
        for (uint_t d = uint_t(0); d < uint_t(3); ++d) {
            state.center[d] = lower.center[d] + alpha * (upper->center[d] - lower.center[d]);
            state.velocity[d] = (upper->center[d] - lower.center[d]) / dt;
        }
        return state;
    };
    // Reason for edit end: evaluate the optional Phase 3 kinematic sphere from either a constant velocity or piecewise-linear keyframes so the PSM path can update the particle state every timestep.

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
    const uint_t topDampingHeight = boundariesConfig.getParameter<uint_t>("topDampingHeight", uint_t(0));
    const real_t topDampingMaxStrength = boundariesConfig.getParameter<real_t>("topDampingMaxStrength", real_t(0));
    const bool topDampingHorizontal = boundariesConfig.getParameter<bool>("topDampingHorizontal", true);
    const bool topDampingVertical = boundariesConfig.getParameter<bool>("topDampingVertical", true);
    // Reason for edit start: read and validate the new top-only zero-normal-pressure-gradient option so it can only run with setup=Open and the existing slip-like top naming.
    const bool topZeroNormalPressureGradient =
            boundariesConfig.getParameter<bool>("topZeroNormalPressureGradient", false);
    const std::string topBoundaryTypeName =
            normalizeConfigToken(boundariesConfig.getParameter<std::string>("topBoundaryType", "symmetry"));
    const bool topBoundaryIsSlipLike =
            topBoundaryTypeName == "symmetry" || topBoundaryTypeName == "slip" || topBoundaryTypeName == "freeslip";
    if (topZeroNormalPressureGradient) {
        WALBERLA_CHECK(boundarySetup.environmentSetup() == EnvironmentSetup::Open,
                       "Boundaries::topZeroNormalPressureGradient requires Boundaries::setup = Open.");
        WALBERLA_CHECK(topBoundaryIsSlipLike,
                       "Boundaries::topZeroNormalPressureGradient requires Boundaries::topBoundaryType to be slip, symmetry, or freeslip.");
    }
    // Reason for edit end: read and validate the new top-only zero-normal-pressure-gradient option so it can only run with setup=Open and the existing slip-like top naming.

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
    const BlockDataID meanVelocityOutputFieldCpuID = walberla::field::addToStorage<VectorField_T>(blocks, "mean velocity output field CPU", real_t(0), layout, fieldGhostLayers, allocator);
    const BlockDataID meanVelocityWfbFieldCpuID = walberla::field::addToStorage<VectorField_T>(blocks, "mean velocity wfb field CPU", real_t(0), layout, fieldGhostLayers, allocator);
    const BlockDataID meanVelocityDampingFieldCpuID = walberla::field::addToStorage<VectorField_T>(blocks, "mean velocity damping field CPU", real_t(0), layout, fieldGhostLayers, allocator);
    const BlockDataID pdfFieldCpuID = walberla::lbm_generated::addPdfFieldToStorage(
            blocks, "pdf field CPU", storageSpec, fieldGhostLayers, layout,
            walberla::Set<walberla::SUID>::emptySet(), walberla::Set<walberla::SUID>::emptySet(), allocator);
    const BlockDataID sumOfSquaresFieldCpuID = walberla::field::addToStorage<SecondOrderTensorField_T>(blocks, "sum of squares field CPU", real_t(0), layout, fieldGhostLayers, allocator);

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
    const uint_t welfordInterval = walberlaEnv.config()->getOneBlock("Output").getParameter<uint_t>("welfordInterval", uint_t(0));
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
    initialiser->setViaVelocityField<VectorField_T, PdfSetter_T>(meanVelocityOutputFieldCpuID, setter);
    initialiser->setViaVelocityField<VectorField_T, PdfSetter_T>(meanVelocityWfbFieldCpuID, setter);
    initialiser->setViaVelocityField<VectorField_T, PdfSetter_T>(meanVelocityDampingFieldCpuID, setter);

    WALBERLA_LOG_INFO_ON_ROOT("Add GPU fields...")

    const BlockDataID densityFieldGpuID = walberla::gpu::addGPUFieldToStorage<ScalarField_T>(blocks, densityFieldCpuID, "density field GPU", true);
    const BlockDataID eddyViscosityFieldGpuID = walberla::gpu::addGPUFieldToStorage<ScalarField_T>(blocks, eddyViscosityFieldCpuID, "eddy viscosity field GPU", true);
    const BlockDataID forceFieldGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T>(blocks, forceFieldCpuID, "force field GPU", true);
    const BlockDataID meanVelocityOutputFieldGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T>(blocks, meanVelocityOutputFieldCpuID, "mean velocity output field GPU", true);
    const BlockDataID meanVelocityWfbFieldGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T>(blocks, meanVelocityWfbFieldCpuID, "mean velocity wfb field GPU", true);
    const BlockDataID meanVelocityDampingFieldGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T>(blocks, meanVelocityDampingFieldCpuID, "mean velocity damping field GPU", true);
    const BlockDataID omegaFieldGpuID = walberla::gpu::addGPUFieldToStorage<ScalarField_T>(blocks, omegaFieldCpuID, "omega field GPU", true);
    const BlockDataID pdfFieldGpuID = walberla::lbm_generated::addGPUPdfFieldToStorage<PdfField_T, StorageSpecification_T>(blocks, pdfFieldCpuID, storageSpec, "pdf field GPU");
    const BlockDataID sumOfSquaresFieldGpuID = walberla::gpu::addGPUFieldToStorage<SecondOrderTensorField_T>(blocks, sumOfSquaresFieldCpuID, "sum of squares field GPU", true);
    const BlockDataID velocityFieldGpuID = walberla::gpu::addGPUFieldToStorage<VectorField_T>(blocks, velocityFieldCpuID, "velocity field GPU", true);
    // Reason for edit start: collect the GPU-backed scalar and vector fields that can be sampled by prm-defined line outputs without needing extra CPU copies.
    std::map<output::Fields::Types, BlockDataID> lineOutputFieldMap{};
    lineOutputFieldMap[output::Fields::DENSITY] = densityFieldGpuID;
    lineOutputFieldMap[output::Fields::VELOCITY] = velocityFieldGpuID;
    lineOutputFieldMap[output::Fields::MEAN_VELOCITY_WFB] = meanVelocityWfbFieldGpuID;
    lineOutputFieldMap[output::Fields::MEAN_VELOCITY_OUTPUT] = meanVelocityOutputFieldGpuID;
    lineOutputFieldMap[output::Fields::FORCE] = forceFieldGpuID;
    lineOutputFieldMap[output::Fields::OMEGA] = omegaFieldGpuID;
    lineOutputFieldMap[output::Fields::EDDY_VISCOSITY] = eddyViscosityFieldGpuID;
    // Reason for edit end: collect the GPU-backed scalar and vector fields that can be sampled by prm-defined line outputs without needing extra CPU copies.

    SweepCollection_T sweepCollection(blocks, densityFieldGpuID, eddyViscosityFieldGpuID, forceFieldGpuID, omegaFieldGpuID,
                                      pdfFieldGpuID, velocityFieldGpuID, omega);

    // Reason for edit start: allocate the optional Phase 3 MesaPD and PSM runtime objects beside the existing GPU solver fields so one kinematic sphere can move through the unchanged urban ABL path.
    using ParticleAccessor_T = walberla::mesa_pd::data::ParticleAccessorWithShape;
    using PSMParticleAndVolumeFractionSoA_T =
            walberla::lbm_mesapd_coupling::psm::gpu::ParticleAndVolumeFractionSoA_T<1>;
    using PSMSweepCollection_T = walberla::lbm_mesapd_coupling::psm::gpu::PSMSweepCollection<
            ParticleAccessor_T,
            walberla::lbm_mesapd_coupling::GlobalParticlesSelector,
            1>;

    std::shared_ptr<walberla::mesa_pd::data::ParticleStorage> psmParticleStorage{};
    std::shared_ptr<walberla::mesa_pd::data::ShapeStorage> psmShapeStorage{};
    std::shared_ptr<ParticleAccessor_T> psmParticleAccessor{};
    std::unique_ptr<PSMParticleAndVolumeFractionSoA_T> psmParticleAndVolumeFractionSoA{};
    std::unique_ptr<PSMSweepCollection_T> psmSweepCollection{};
    std::unique_ptr<walberla::pystencils::waLBerlaABLPSM_Sweep> psmSweep{};
    std::unique_ptr<walberla::pystencils::waLBerlaABLPSM_InitializeDomainForPSM> psmPdfInitializer{};
    walberla::lbm_mesapd_coupling::GlobalParticlesSelector psmGlobalParticleSelector{};
    walberla::id_t psmSphereUID(0);

    if (psmCouplingEnabled) {
        const walberla::AABB & domainAABB = blocks->getDomain();
        const auto validateSphereCenter = [&](const walberla::Vector3<real_t> & center, const char * description) {
            for (uint_t d = uint_t(0); d < uint_t(3); ++d) {
                WALBERLA_CHECK(center[d] - movingBodyRuntimeConfig.body.radius >= domainAABB.min(d) &&
                               center[d] + movingBodyRuntimeConfig.body.radius <= domainAABB.max(d),
                               description)
            }
        };

        if (movingBodyRuntimeConfig.body.trajectoryPoints.empty()) {
            validateSphereCenter(movingBodyRuntimeConfig.body.initialPosition,
                                 "Phase 3 moving sphere must start fully inside the simulation domain.");
        } else {
            for (const auto & point : movingBodyRuntimeConfig.body.trajectoryPoints) {
                validateSphereCenter(point.center,
                                     "Phase 3 moving sphere trajectory points must stay fully inside the simulation domain.");
            }
        }

        psmParticleStorage = std::make_shared<walberla::mesa_pd::data::ParticleStorage>(1);
        psmShapeStorage = std::make_shared<walberla::mesa_pd::data::ShapeStorage>();
        psmParticleAccessor = std::make_shared<ParticleAccessor_T>(psmParticleStorage, psmShapeStorage);

        const auto sphereShape =
                psmShapeStorage->create<walberla::mesa_pd::data::Sphere>(movingBodyRuntimeConfig.body.radius);
        const MovingSphereState initialMovingSphereState = evaluateMovingSphereState(uint_t(0));
        const walberla::mesa_pd::Vec3 movingSpherePosition(
                initialMovingSphereState.center[0],
                initialMovingSphereState.center[1],
                initialMovingSphereState.center[2]);
        const walberla::mesa_pd::Vec3 movingSphereVelocity(
                initialMovingSphereState.velocity[0],
                initialMovingSphereState.velocity[1],
                initialMovingSphereState.velocity[2]);
        walberla::mesa_pd::data::ParticleStorage::Particle&& movingSphere =
                *psmParticleStorage->create(true);
        movingSphere.setPosition(movingSpherePosition);
        movingSphere.setLinearVelocity(movingSphereVelocity);
        movingSphere.setInteractionRadius(movingBodyRuntimeConfig.body.radius);
        movingSphere.setOwner(walberla::mpi::MPIManager::instance()->rank());
        movingSphere.setShapeID(sphereShape);
        psmSphereUID = movingSphere.getUid();

        psmParticleAndVolumeFractionSoA =
                std::make_unique<PSMParticleAndVolumeFractionSoA_T>(blocks, omega);
        psmSweepCollection = std::make_unique<PSMSweepCollection_T>(
                blocks,
                psmParticleAccessor,
                psmGlobalParticleSelector,
                *psmParticleAndVolumeFractionSoA,
                walberla::Vector3<uint_t>(uint_t(8), uint_t(8), uint_t(8)));
        psmPdfInitializer =
                std::make_unique<walberla::pystencils::waLBerlaABLPSM_InitializeDomainForPSM>(
                        psmParticleAndVolumeFractionSoA->BsFieldID,
                        psmParticleAndVolumeFractionSoA->BFieldID,
                        densityFieldGpuID,
                        forceFieldGpuID,
                        psmParticleAndVolumeFractionSoA->particleVelocitiesFieldID,
                        pdfFieldGpuID,
                        velocityFieldGpuID);
        psmSweep = std::make_unique<walberla::pystencils::waLBerlaABLPSM_Sweep>(
                psmParticleAndVolumeFractionSoA->BsFieldID,
                psmParticleAndVolumeFractionSoA->BFieldID,
                densityFieldGpuID,
                eddyViscosityFieldGpuID,
                forceFieldGpuID,
                omegaFieldGpuID,
                psmParticleAndVolumeFractionSoA->particleForcesFieldID,
                psmParticleAndVolumeFractionSoA->particleVelocitiesFieldID,
                pdfFieldGpuID,
                velocityFieldGpuID,
                omega);

        WALBERLA_LOG_INFO_ON_ROOT("Phase 3 kinematic sphere initial state: center = "
                                  << initialMovingSphereState.center
                                  << ", velocity = " << initialMovingSphereState.velocity
                                  << ", radius = " << movingBodyRuntimeConfig.body.radius)
    }
    // Reason for edit end: allocate the optional Phase 3 MesaPD and PSM runtime objects beside the existing GPU solver fields so one kinematic sphere can move through the unchanged urban ABL path.

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

    walberla::pystencils::waLBerlaABL_WelfordWFB welfordTopDampingSweep(
            meanVelocityDampingFieldGpuID, velocityFieldGpuID,
            real_t(0));
    auto welfordTopDampingLambda = [&welfordTopDampingSweep](walberla::IBlock * block) { welfordTopDampingSweep(block); };

        walberla::pystencils::waLBerlaABL_SoSResetter welfordOutputSosResetter(sumOfSquaresFieldGpuID);
        walberla::pystencils::waLBerlaABL_WelfordOutput welfordOutputSweep(
            meanVelocityOutputFieldGpuID, sumOfSquaresFieldGpuID, velocityFieldGpuID,
            real_t(0));
        auto welfordOutputLambda = [&welfordOutputSweep](walberla::IBlock * block) { welfordOutputSweep(block); };

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
    const real_t vtkSamplingResolution = vtkConfig.getParameter<real_t>("samplingResolution", real_t(-1));
    // Reason for edit start: store shared dynamic AABB filter state and version-1 piecewise-linear trajectory keyframes so copied VTK filters all observe the same runtime-updated bounds.
    struct DynamicAABBKeyframe {
        uint_t step;
        walberla::Vector3<real_t> center;
        walberla::Vector3<real_t> size;
    };
    struct DynamicAABBMotion {
        std::shared_ptr<output::DynamicAABBInclusionFilter::State> state;
        std::vector<DynamicAABBKeyframe> keyframes;
        bool clampToDomain;
    };
    std::vector<DynamicAABBMotion> dynamicAABBMotions{};
    // Reason for edit end: store shared dynamic AABB filter state and version-1 piecewise-linear trajectory keyframes so copied VTK filters all observe the same runtime-updated bounds.
    std::shared_ptr<walberla::vtk::VTKOutput> fieldVTKOutput{nullptr};
    if(vtkWriteFrequency > 0) {
        fieldVTKOutput = walberla::vtk::createVTKOutput_BlockData(
                *blocks, //structured block storage
                "abl_gpu_uniform", // Identifier
            uint_t(1), //writefrequency
                vtkGhostLayers,
                false, //forcePVTU
                vtkBaseFolder,
                vtkExecutionFolder,
                true, //continuousNumbering
                true, // binaryOutput
                true, // littleEndianOutput
                true, // useMPIIO
                uint_t(0), //Initial Execution Count
                false, //amrFileFormat
                false); //oneFilePerProcess
        fieldVTKOutput->setSamplingResolution(vtkSamplingResolution);
        if(vtkConfig.isDefined("samplingDx")) {
            const real_t vtkSamplingDx = vtkConfig.getParameter<real_t>("samplingDx", real_t(-1));
            const real_t vtkSamplingDy = vtkConfig.getParameter<real_t>("samplingDy", real_t(-1));
            const real_t vtkSamplingDz = vtkConfig.getParameter<real_t>("samplingDz", real_t(-1));
            fieldVTKOutput->setSamplingResolution(vtkSamplingDx, vtkSamplingDy, vtkSamplingDz);
        }
        // fieldVTKOutput->addCellDataWriter(std::make_shared<walberla::field::VTKWriter<ScalarField_T>>(densityFieldCpuID, "Density"));
        fieldVTKOutput->addCellDataWriter(std::make_shared<walberla::field::VTKWriter<VectorField_T>>(velocityFieldCpuID, "Velocity"));
        fieldVTKOutput->addCellDataWriter(std::make_shared<walberla::field::VTKWriter<VectorField_T>>(meanVelocityOutputFieldCpuID, "MeanVelocityOutput"));
        fieldVTKOutput->addCellDataWriter(std::make_shared<walberla::field::VTKWriter<SecondOrderTensorField_T>>(sumOfSquaresFieldCpuID, "SumOfSquares"));
        // fieldVTKOutput->addCellDataWriter(std::make_shared<walberla::field::VTKWriter<VectorField_T>>(forceFieldCpuID, "Force"));
        // fieldVTKOutput->addCellDataWriter(std::make_shared<walberla::field::VTKWriter<ScalarField_T>>(eddyViscosityFieldCpuID, "EddyViscosity"));
        // fieldVTKOutput->addCellDataWriter(std::make_shared<walberla::field::VTKWriter<ScalarField_T>>(omegaFieldCpuID, "Omega"));
        // fieldVTKOutput->addCellDataWriter(std::make_shared<walberla::field::VTKWriter<FlagField_T>>(flagFieldID, "Flag"));
        // Reason for edit start: parse GeneratedApp-style inclusion filters from the prm file and attach them to the existing hardcoded writer so users can switch between whole-domain and filtered VTK output without losing the current default behavior.
        auto addDomainFilter = [&](walberla::vtk::ChainedFilter * chainedFilter = nullptr) {
            walberla::field::FlagFieldCellFilter<FlagField_T> fluidFilter(flagFieldID);
            fluidFilter.addFlag(FluidFlagUID);
            if(chainedFilter != nullptr) {
                chainedFilter->addFilter(fluidFilter);
            } else {
                fieldVTKOutput->addCellInclusionFilter(fluidFilter);
            }
        };
        auto addAABBFilters = [&](const walberla::Config::BlockHandle & filterConfig,
                                  walberla::vtk::ChainedFilter * chainedFilter = nullptr) {
            std::vector<walberla::Config::BlockHandle> aabbFilterBlocks;
            filterConfig.getBlocks("AABB", aabbFilterBlocks);
            for(const auto & filter : aabbFilterBlocks) {
                const walberla::Vector3<real_t> min = filter.getParameter<walberla::Vector3<real_t>>("min");
                const walberla::Vector3<real_t> max = filter.getParameter<walberla::Vector3<real_t>>("max");
                // auto min_local = min;
                // auto max_local = max;
                // min_local[0] += 0;  //test
                // max_local[0] += 0;   //test
                walberla::vtk::AABBCellFilter aabbFilter({min, max}); 
                // walberla::vtk::AABBCellFilter aabbFilter({min_local, max_local}); // commented for trying dynamic min max
                if(chainedFilter != nullptr) {
                    chainedFilter->addFilter(aabbFilter);
                } else {
                    fieldVTKOutput->addCellInclusionFilter(aabbFilter);
                }
            }
        };
        // Reason for edit start: parse a version-1 DynamicAABB prm block into shared-state filters driven by piecewise-linear center/size keyframes with domain clamping.
        auto addDynamicAABBFilters = [&](const walberla::Config::BlockHandle & filterConfig,
                                         walberla::vtk::ChainedFilter * chainedFilter = nullptr) {
            std::vector<walberla::Config::BlockHandle> dynamicAABBFilterBlocks;
            filterConfig.getBlocks("DynamicAABB", dynamicAABBFilterBlocks);
            for(const auto & filter : dynamicAABBFilterBlocks) {
                const std::string trajectoryType =
                        normalizeConfigToken(filter.getParameter<std::string>("trajectoryType", "piecewise_linear"));
                const std::string limitMode =
                        normalizeConfigToken(filter.getParameter<std::string>("limitMode", "clamp"));
                WALBERLA_CHECK(trajectoryType == "piecewise_linear",
                               "DynamicAABB currently supports only trajectoryType = piecewise_linear.");
                WALBERLA_CHECK(limitMode == "clamp",
                               "DynamicAABB currently supports only limitMode = clamp.");

                const walberla::Vector3<real_t> center0 =
                        filter.getParameter<walberla::Vector3<real_t>>("center0");
                const walberla::Vector3<real_t> size0 =
                        filter.getParameter<walberla::Vector3<real_t>>("size0");

                for(uint_t d = 0; d < uint_t(3); ++d) {
                    WALBERLA_CHECK_GREATER(size0[d], real_t(0),
                                           "DynamicAABB size0 components must be positive.");
                }

                auto min0 = center0;
                auto max0 = center0;
                for(uint_t d = 0; d < uint_t(3); ++d) {
                    const real_t halfSize = real_t(0.5) * size0[d];
                    min0[d] -= halfSize;
                    max0[d] += halfSize;
                }

                auto state = std::make_shared<output::DynamicAABBInclusionFilter::State>(min0, max0);

                DynamicAABBMotion motion{ state, {}, true };
                motion.keyframes.push_back(DynamicAABBKeyframe{ uint_t(0), center0, size0 });

                std::vector<walberla::Config::BlockHandle> pointBlocks;
                filter.getBlocks("point", pointBlocks);
                std::sort(pointBlocks.begin(), pointBlocks.end(),
                          [](const walberla::Config::BlockHandle & lhs, const walberla::Config::BlockHandle & rhs) {
                              return lhs.getParameter<uint_t>("step") < rhs.getParameter<uint_t>("step");
                          });

                auto currentSize = size0;
                for(const auto & point : pointBlocks) {
                    const uint_t step = point.getParameter<uint_t>("step");
                    const walberla::Vector3<real_t> center =
                            point.getParameter<walberla::Vector3<real_t>>("center");
                    if(point.isDefined("size")) {
                        currentSize = point.getParameter<walberla::Vector3<real_t>>("size");
                        for(uint_t d = 0; d < uint_t(3); ++d) {
                            WALBERLA_CHECK_GREATER(currentSize[d], real_t(0),
                                                   "DynamicAABB point size components must be positive.");
                        }
                    }

                    DynamicAABBKeyframe keyframe{ step, center, currentSize };
                    if(step == motion.keyframes.back().step) {
                        motion.keyframes.back() = keyframe;
                    } else {
                        motion.keyframes.push_back(keyframe);
                    }
                }

                dynamicAABBMotions.push_back(motion);

                output::DynamicAABBInclusionFilter dynamicAABBFilter(state);
                if(chainedFilter != nullptr) {
                    chainedFilter->addFilter(dynamicAABBFilter);
                } else {
                    fieldVTKOutput->addCellInclusionFilter(dynamicAABBFilter);
                }
            }
        };
        // Reason for edit end: parse a version-1 DynamicAABB prm block into shared-state filters driven by piecewise-linear center/size keyframes with domain clamping.
        auto addPlaneFilters = [&](const walberla::Config::BlockHandle & filterConfig,
                                   walberla::vtk::ChainedFilter * chainedFilter = nullptr) {
            std::vector<walberla::Config::BlockHandle> planeFilterBlocks;
            filterConfig.getBlocks("Plane", planeFilterBlocks);
            for(const auto & filter : planeFilterBlocks) {
                const walberla::Vector3<real_t> point = filter.getParameter<walberla::Vector3<real_t>>("point");
                const walberla::Vector3<real_t> normal = filter.getParameter<walberla::Vector3<real_t>>("normal");
                const real_t maxDistance = filter.getParameter<real_t>("maxDistance", real_t(0.5));
                output::PlaneInclusionFilter planeFilter(point, normal, maxDistance);
                if(chainedFilter != nullptr) {
                    chainedFilter->addFilter(planeFilter);
                } else {
                    fieldVTKOutput->addCellInclusionFilter(planeFilter);
                }
            }
        };
        if(vtkConfig.getNumBlocks("inclusion_filters")) {
            WALBERLA_LOG_INFO_ON_ROOT("ABLApp_gpu_uniform VTK inclusion filters enabled")
            auto inclusionBlock = vtkConfig.getOneBlock("inclusion_filters");
            if(inclusionBlock.isDefined("DomainFilter")) {
                addDomainFilter();
            }
            if(inclusionBlock.getNumBlocks("AABB")) {
                addAABBFilters(inclusionBlock);
            }
            if(inclusionBlock.getNumBlocks("DynamicAABB")) {
                addDynamicAABBFilters(inclusionBlock);
            }
            if(inclusionBlock.getNumBlocks("Plane")) {
                addPlaneFilters(inclusionBlock);
            }
            std::vector<walberla::Config::BlockHandle> combineFilterBlocks;
            inclusionBlock.getBlocks("combine", combineFilterBlocks);
            for(const auto & combineFilter : combineFilterBlocks) {
                walberla::vtk::ChainedFilter chainedFilter{};
                if(combineFilter.isDefined("DomainFilter")) {
                    addDomainFilter(&chainedFilter);
                }
                if(combineFilter.getNumBlocks("AABB")) {
                    addAABBFilters(combineFilter, &chainedFilter);
                }
                if(combineFilter.getNumBlocks("DynamicAABB")) {
                    addDynamicAABBFilters(combineFilter, &chainedFilter);
                }
                if(combineFilter.getNumBlocks("Plane")) {
                    addPlaneFilters(combineFilter, &chainedFilter);
                }
                fieldVTKOutput->addCellInclusionFilter(chainedFilter);
            }
        }
        // Reason for edit end: parse GeneratedApp-style inclusion filters from the prm file and attach them to the existing hardcoded writer so users can switch between whole-domain and filtered VTK output without losing the current default behavior.
    }

    // Reason for edit start: evaluate version-1 DynamicAABB piecewise-linear trajectories each timestep and clamp the resulting window to the simulation domain.
    auto updateDynamicAABBFilters = [&]() {
        if(dynamicAABBMotions.empty()) {
            return;
        }

        static uint_t dynamicAABBStepCounter = uint_t(0);
        const walberla::AABB & domainAABB = blocks->getDomain();

        for(auto & motion : dynamicAABBMotions) {
            WALBERLA_CHECK(!motion.keyframes.empty(), "DynamicAABB requires at least one keyframe.");

            walberla::Vector3<real_t> center = motion.keyframes.back().center;
            walberla::Vector3<real_t> size = motion.keyframes.back().size;

            const auto upper = std::lower_bound(
                    motion.keyframes.begin(), motion.keyframes.end(), dynamicAABBStepCounter,
                    [](const DynamicAABBKeyframe & keyframe, const uint_t step) { return keyframe.step < step; });

            if(upper == motion.keyframes.begin()) {
                center = upper->center;
                size = upper->size;
            } else if(upper != motion.keyframes.end()) {
                if(upper->step == dynamicAABBStepCounter) {
                    center = upper->center;
                    size = upper->size;
                } else {
                    const auto & lower = *(upper - 1);
                    const real_t intervalLength = real_t(upper->step - lower.step);
                    const real_t alpha =
                            intervalLength > real_t(0)
                            ? real_t(dynamicAABBStepCounter - lower.step) / intervalLength
                            : real_t(0);
                    for(uint_t d = 0; d < uint_t(3); ++d) {
                        center[d] = lower.center[d] + alpha * (upper->center[d] - lower.center[d]);
                        size[d] = lower.size[d] + alpha * (upper->size[d] - lower.size[d]);
                    }
                }
            }

            for(uint_t d = 0; d < uint_t(3); ++d) {
                WALBERLA_CHECK_GREATER(size[d], real_t(0),
                                       "DynamicAABB interpolated size components must stay positive.");
                const real_t domainExtent = domainAABB.max(d) - domainAABB.min(d);
                if(motion.clampToDomain) {
                    size[d] = std::min(size[d], domainExtent);
                    const real_t halfSize = real_t(0.5) * size[d];
                    center[d] = std::clamp(center[d], domainAABB.min(d) + halfSize, domainAABB.max(d) - halfSize);
                }
            }

            auto min = center;
            auto max = center;
            for(uint_t d = 0; d < uint_t(3); ++d) {
                const real_t halfSize = real_t(0.5) * size[d];
                min[d] -= halfSize;
                max[d] += halfSize;
            }
            motion.state->setBounds(min, max);
        }

        ++dynamicAABBStepCounter;
    };
    // Reason for edit end: evaluate version-1 DynamicAABB piecewise-linear trajectories each timestep and clamp the resulting window to the simulation domain.

    auto writeVTK = [&]() {
        if(!fieldVTKOutput) {
            return;
        }

        static uint_t vtkStepCounter = uint_t(0);
        if(vtkStepCounter < vtkStartTimestep) {
            ++vtkStepCounter;
            return;
        }

        if(vtkWriteFrequency == uint_t(0)) {
            ++vtkStepCounter;
            return;
        }

        // Manual cadence control: only perform GPU->CPU copies and synchronization on real output steps.
        const uint_t stepsSinceStart = vtkStepCounter - vtkStartTimestep;
        if(stepsSinceStart % vtkWriteFrequency != uint_t(0)) {
            ++vtkStepCounter;
            return;
        }

        //walberla::gpu::fieldCpy<ScalarField_T, GPUField_T<real_t>>(blocks, densityFieldCpuID, densityFieldGpuID);
        walberla::gpu::fieldCpy<VectorField_T, GPUField_T<real_t>>(blocks, velocityFieldCpuID, velocityFieldGpuID);
        //walberla::gpu::fieldCpy<VectorField_T, GPUField_T<real_t>>(blocks, meanVelocityOutputFieldCpuID, meanVelocityOutputFieldGpuID);
        //walberla::gpu::fieldCpy<SecondOrderTensorField_T, GPUField_T<real_t>>(blocks, sumOfSquaresFieldCpuID, sumOfSquaresFieldGpuID);
        //walberla::gpu::fieldCpy<VectorField_T, GPUField_T<real_t>>(blocks, forceFieldCpuID, forceFieldGpuID);
        //walberla::gpu::fieldCpy<ScalarField_T, GPUField_T<real_t>>(blocks, eddyViscosityFieldCpuID, eddyViscosityFieldGpuID);
        //walberla::gpu::fieldCpy<ScalarField_T, GPUField_T<real_t>>(blocks, omegaFieldCpuID, omegaFieldGpuID);
        // Keep synchronization for correctness: host VTK writer must see completed GPU->CPU copies.
        cudaDeviceSynchronize();

        fieldVTKOutput->write();
        ++vtkStepCounter;
    };
    // Reason for edit start: instantiate the optional line-output helper so prm-defined probe lines can run beside or instead of VTK output without changing the existing hardcoded writer flow.
    std::shared_ptr<LineOutput_T> lineOutput{nullptr};
    if(outputConfig.getNumBlocks("LineOutput")) {
        WALBERLA_LOG_INFO_ON_ROOT("ABLApp_gpu_uniform line output enabled")
        lineOutput = std::make_shared<LineOutput_T>(outputConfig, blocks, &timeloop, lineOutputFieldMap);
    }
    auto writeLineOutput = [&]() {
        if(!lineOutput) {
            return;
        }
        lineOutput->write();
    };
    // Reason for edit end: instantiate the optional line-output helper so prm-defined probe lines can run beside or instead of VTK output without changing the existing hardcoded writer flow.
// VTK output (GPU-safe): copy selected fields to CPU before writing
    walberla::wind::FlowDriverCollection flowDriver(blocks, &timeloop, globalConfig, domainSetup, forceFieldGpuID, velocityFieldGpuID, fieldGhostLayers);
    for(auto & block : *blocks) { flowDriver(&block); }
    // Reason for edit start: construct the top-only zero-normal-pressure-gradient sweep as an opt-in app-level extension that leaves the bottom WFB logic unchanged.
    boundary::TopSlipZeroGradientBoundary<PdfGPUField_T, GPUField_T<real_t>> topSlipZeroGradientBoundary(
            blocks,
            pdfFieldGpuID,
            forceFieldGpuID,
            domainSetup.domainSize_[2],
            topZeroNormalPressureGradient);
    if(topSlipZeroGradientBoundary.isEnabled()) {
        WALBERLA_LOG_INFO_ON_ROOT("Top slip zero-gradient pressure handling enabled")
    }
    // Reason for edit end: construct the top-only zero-normal-pressure-gradient sweep as an opt-in app-level extension that leaves the bottom WFB logic unchanged.
    // Reason for edit start: replace the generated top FreeSlip handling with a custom boundary sweep when the top zero-normal-pressure-gradient option is active, so the top boundary is not applied twice while all other boundary operators keep their original ordering.
    std::function<void(walberla::IBlock *)> boundaryHandlingSweep = boundaryCollection.getSweep();
    if(topSlipZeroGradientBoundary.isEnabled()) {
        boundaryHandlingSweep = [&boundaryCollection](walberla::IBlock * block) {
            boundaryCollection.waLBerlaABL_NoSlipObject->run(block);
            boundaryCollection.waLBerlaABL_WFBObject->run(block);
            boundaryCollection.waLBerlaABL_UniformUBBObject->run(block);
            boundaryCollection.waLBerlaABL_LogLawUBBObject->run(block);
            boundaryCollection.waLBerlaABL_OutflowObject->run(block);
            boundaryCollection.waLBerlaABL_TopOutflowObject->run(block);
        };
        WALBERLA_LOG_INFO_ON_ROOT("Generated top FreeSlip handling disabled while top slip zero-gradient pressure handling is active")
    }
    // Reason for edit end: replace the generated top FreeSlip handling with a custom boundary sweep when the top zero-normal-pressure-gradient option is active, so the top boundary is not applied twice while all other boundary operators keep their original ordering.

    damping::TopDampingZone topDamping(
            blocks,
            forceFieldGpuID,
            velocityFieldGpuID,
            meanVelocityDampingFieldGpuID,
            domainSetup.domainSize_[2],
            topDampingHeight,
            topDampingMaxStrength,
            topDampingHorizontal,
            topDampingVertical);
    if(topDamping.isEnabled()) {
        WALBERLA_LOG_INFO_ON_ROOT("Top damping enabled: height = " << topDampingHeight
                                  << ", maxStrength = " << topDampingMaxStrength
                                  << ", dampHorizontal = " << topDampingHorizontal
                                  << ", dampVertical = " << topDampingVertical)
        for(auto & block : *blocks) { topDamping(&block); }
    }

    // Reason for edit start: update the optional Phase 3 kinematic sphere from the prm trajectory before mapping so MesaPD and the PSM sweeps always see the same current particle state.
    const auto applyMovingSphereState = [&](const uint_t step, const bool logState = false) {
        if(!psmCouplingEnabled) {
            return;
        }

        const MovingSphereState state = evaluateMovingSphereState(step);
        const walberla::AABB & domainAABB = blocks->getDomain();
        for(uint_t d = uint_t(0); d < uint_t(3); ++d) {
            WALBERLA_CHECK(state.center[d] - movingBodyRuntimeConfig.body.radius >= domainAABB.min(d) &&
                           state.center[d] + movingBodyRuntimeConfig.body.radius <= domainAABB.max(d),
                           "Phase 3 moving sphere left the simulation domain.")
        }

        const size_t sphereIdx = psmParticleAccessor->uidToIdx(psmSphereUID);
        WALBERLA_CHECK(sphereIdx != psmParticleAccessor->getInvalidIdx(),
                       "Phase 3 moving sphere particle could not be found in the MesaPD accessor.")

        psmParticleAccessor->setPosition(
                sphereIdx,
                walberla::mesa_pd::Vec3(state.center[0], state.center[1], state.center[2]));
        psmParticleAccessor->setLinearVelocity(
                sphereIdx,
                walberla::mesa_pd::Vec3(state.velocity[0], state.velocity[1], state.velocity[2]));

        if(logState) {
            WALBERLA_LOG_INFO_ON_ROOT("Phase 3 moving sphere applied at step " << step
                                      << ": center = " << state.center
                                      << ", velocity = " << state.velocity)
        }
    };
    // Reason for edit end: update the optional Phase 3 kinematic sphere from the prm trajectory before mapping so MesaPD and the PSM sweeps always see the same current particle state.

    // Reason for edit start: map the Phase 3 kinematic sphere and reinitialize the surrounding PDFs once before the first timestep so the optional PSM path starts from a consistent fluid state.
    if(psmCouplingEnabled) {
        applyMovingSphereState(uint_t(0), true);
        for(auto & block : *blocks) {
            psmSweepCollection->particleMappingSweep(&block);
        }
        for(auto & block : *blocks) {
            psmSweepCollection->setParticleVelocitiesSweep(&block);
            (*psmPdfInitializer)(&block);
        }
        cudaDeviceSynchronize();
    }
    // Reason for edit end: map the Phase 3 kinematic sphere and reinitialize the surrounding PDFs once before the first timestep so the optional PSM path starts from a consistent fluid state.

    timeloop.add() << walberla::BeforeFunction(communication->getCommunicateFunctor(), "Field communication")
                   << walberla::BeforeFunction([&boundarySetup, &shiftedPeriodicity]() {
                       if(boundarySetup.inflowType() == InflowSetup::ShiftedPeriodic) shiftedPeriodicity();
                   }, "Shifted periodicity")
                   << walberla::Sweep(boundaryHandlingSweep, "Boundary handling");
    // Reason for edit start: run the top-only zero-normal-pressure-gradient sweep immediately after the boundary handling sweep so it becomes the only active top treatment when the generated FreeSlip path has been skipped.
    if(topSlipZeroGradientBoundary.isEnabled()) {
        timeloop.add() << walberla::Sweep(topSlipZeroGradientBoundary.getSweep(), "Top slip zero-gradient pressure handling");
    }
    // Reason for edit end: run the top-only zero-normal-pressure-gradient sweep immediately after the boundary handling sweep so it becomes the only active top treatment when the generated FreeSlip path has been skipped.

    // Reason for edit start: keep the legacy stream-collide path untouched when PSM is off and swap in the synchronized Phase 3 kinematic-sphere sweep sequence only when the optional coupling is enabled.
    if(psmCouplingEnabled) {
        const auto updateMovingSphereTrajectory = [&]() {
            applyMovingSphereState(timeloop.getCurrentTimeStep());
        };
        const std::function<void(walberla::IBlock *)> psmParticleMappingSweep =
                [&psmSweepCollection](walberla::IBlock * block) { psmSweepCollection->particleMappingSweep(block); };
        const std::function<void(walberla::IBlock *)> psmSetParticleVelocitiesSweep =
                [&psmSweepCollection](walberla::IBlock * block) { psmSweepCollection->setParticleVelocitiesSweep(block); };
        const std::function<void(walberla::IBlock *)> psmStreamCollideSweep = psmSweep->getSweep();
        const std::function<void(walberla::IBlock *)> psmReduceParticleForcesSweep =
                [&psmSweepCollection](walberla::IBlock * block) { psmSweepCollection->reduceParticleForcesSweep(block); };

        timeloop.add() << walberla::BeforeFunction(updateMovingSphereTrajectory, "Update moving sphere trajectory")
                       << walberla::Sweep(
                                walberla::lbm_mesapd_coupling::psm::gpu::deviceSyncWrapper(psmParticleMappingSweep),
                                "Particle mapping");
        timeloop.add() << walberla::Sweep(
                                walberla::lbm_mesapd_coupling::psm::gpu::deviceSyncWrapper(psmSetParticleVelocitiesSweep),
                                "Set particle velocities");
        timeloop.add() << walberla::Sweep(
                                walberla::lbm_mesapd_coupling::psm::gpu::deviceSyncWrapper(psmStreamCollideSweep),
                                "ABL PSM stream-collide");
        timeloop.add() << walberla::Sweep(
                                walberla::lbm_mesapd_coupling::psm::gpu::deviceSyncWrapper(psmReduceParticleForcesSweep),
                                "Reduce particle forces")
                       << walberla::AfterFunction(updateDynamicAABBFilters, "Dynamic AABB filter update")
                       << walberla::AfterFunction(writeVTK, "VTK output")
                       << walberla::AfterFunction(writeLineOutput, "Line output");
    } else {
        timeloop.add() << walberla::Sweep(sweepCollection.streamCollide(), "LBM stream-collide")
                       << walberla::AfterFunction(updateDynamicAABBFilters, "Dynamic AABB filter update")
                       << walberla::AfterFunction(writeVTK, "VTK output")
                       << walberla::AfterFunction(writeLineOutput, "Line output");
    }
    // Reason for edit end: keep the legacy stream-collide path untouched when PSM is off and swap in the synchronized Phase 3 kinematic-sphere sweep sequence only when the optional coupling is enabled.

    // Temporary benchmark toggle: comment out all Welford timeloop registration
    // blocks instead of deleting them, so we can restore the current version later.

    if(boundarySetup.wallType() == WallSetup::WFB) {
        timeloop.add() << walberla::BeforeFunction([&]() {
                           welfordWFBSweep.setCounter(real_t(welfordWFBSweep.getCounter() + 1));
                       }, "WelfordWFB counter")
                       << walberla::Sweep(welfordWFBLambda, "WelfordWFB sweep");
    }

    timeloop.add() << walberla::BeforeFunction([&]() {
                       if(welfordInterval && (uint_t(welfordOutputSweep.getCounter()) % welfordInterval == 0)) {
                           welfordOutputSweep.setCounter(real_t(1));
                           for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
                               auto dst = blockIt->getData<GPUField_T<real_t>>(meanVelocityOutputFieldGpuID);
                               const auto src = blockIt->getData<GPUField_T<real_t>>(velocityFieldGpuID);
                               gpu::fieldCopy(dst, src);
                               welfordOutputSosResetter(blockIt.get());
                           }
                       } else {
                           welfordOutputSweep.setCounter(real_t(welfordOutputSweep.getCounter() + 1));
                       }
                   }, "WelfordOutput counter")
                   << walberla::Sweep(welfordOutputLambda, "WelfordOutput sweep");

    timeloop.add() << walberla::BeforeFunction([&]() {
                       welfordTopDampingSweep.setCounter(real_t(welfordTopDampingSweep.getCounter() + 1));
                   }, "WelfordTopDamping counter")
                   << walberla::Sweep(welfordTopDampingLambda, "WelfordTopDamping sweep");

    // Temporary benchmark toggle: comment out all Welford timeloop registration
    // blocks instead of deleting them, so we can restore the current version later.
    timeloop.add() << walberla::Sweep(flowDriver, "Setting driving force");
    if(topDamping.isEnabled()) {
        timeloop.add() << walberla::Sweep(topDamping.getSweep(), "Top damping");
    }

    // Reason for edit start: reset the per-timestep hydrodynamic force and torque only when the Phase 3 moving sphere is active so future force evaluation stays clean without changing default ABL runs.
    if(psmCouplingEnabled) {
        timeloop.addFuncAfterTimeStep(
                [psmParticleStorage, psmParticleAccessor]() {
                    psmParticleStorage->forEachParticle(
                            false,
                            walberla::mesa_pd::kernel::SelectAll(),
                            *psmParticleAccessor,
                            walberla::lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel(),
                            *psmParticleAccessor);
                },
                "Reset PSM hydrodynamic force");
    }
    // Reason for edit end: reset the per-timestep hydrodynamic force and torque only when the Phase 3 moving sphere is active so future force evaluation stays clean without changing default ABL runs.

    timeloop.addFuncAfterTimeStep(
            walberla::makeSharedFunctor(
                    walberla::field::makeStabilityChecker<PdfField_T, FlagField_T>(
                            globalConfig, blocks, pdfFieldCpuID, flagFieldID, FluidFlagUID)),
            "LBM stability check");

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
