// Reason for edit start: add a top-only sweep that keeps setup=Open and the bottom WFB path intact while approximating zero normal pressure gradient with slip-style kinematics at the top boundary.
#pragma once

#ifndef TURBINECORE_ABLAPP_GPU_UNIFORM_TOPSLIPZEROGRADIENTBOUNDARY_H
#define TURBINECORE_ABLAPP_GPU_UNIFORM_TOPSLIPZEROGRADIENTBOUNDARY_H

#include <cstddef>
#include <functional>
#include <memory>

#include <blockforest/StructuredBlockForest.h>
#include <core/DataTypes.h>
#include <core/Macros.h>
#include <domain_decomposition/BlockDataID.h>
#include <domain_decomposition/IBlock.h>
#include <gpu/ErrorChecking.h>
#include <gpu/GPUField.h>
#include <gpu/GPUWrapper.h>

namespace turbine_core {
namespace boundary {

namespace internal {

static constexpr int kTopSlipZeroGradientThreadsX = 16;
static constexpr int kTopSlipZeroGradientThreadsY = 16;

template <typename Real_T>
struct GpuFieldView
{
    Real_T * data;
    cell_idx_t xSize;
    cell_idx_t ySize;
    cell_idx_t zSize;
    cell_idx_t fSize;
    cell_idx_t xStride;
    cell_idx_t yStride;
    cell_idx_t zStride;
    cell_idx_t fStride;
    cell_idx_t ghostLayers;

    HOST_DEVICE_PREFIX Real_T & get(cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f) const
    {
        const std::ptrdiff_t idx =
            static_cast<std::ptrdiff_t>(f) * static_cast<std::ptrdiff_t>(fStride) +
            static_cast<std::ptrdiff_t>(x + ghostLayers) * static_cast<std::ptrdiff_t>(xStride) +
            static_cast<std::ptrdiff_t>(y + ghostLayers) * static_cast<std::ptrdiff_t>(yStride) +
            static_cast<std::ptrdiff_t>(z + ghostLayers) * static_cast<std::ptrdiff_t>(zStride);
        return data[idx];
    }
};

template <typename Field_T>
HOST_PREFIX GpuFieldView<typename Field_T::value_type> makeGpuFieldView(Field_T * field)
{
    using Real_T = typename Field_T::value_type;
    const cell_idx_t ghostLayers = field->xOff();

    return GpuFieldView<Real_T>{
        field->dataAt(-ghostLayers, -ghostLayers, -ghostLayers, 0),
        static_cast<cell_idx_t>(field->xSize()),
        static_cast<cell_idx_t>(field->ySize()),
        static_cast<cell_idx_t>(field->zSize()),
        static_cast<cell_idx_t>(field->fSize()),
        static_cast<cell_idx_t>(field->xStride()),
        static_cast<cell_idx_t>(field->yStride()),
        static_cast<cell_idx_t>(field->zStride()),
        static_cast<cell_idx_t>(field->fStride()),
        ghostLayers
    };
}

template <typename Real_T>
HOST_DEVICE_PREFIX void computeMacroscopicState(const GpuFieldView<Real_T> & pdfField,
                                                const GpuFieldView<Real_T> & forceField,
                                                cell_idx_t x,
                                                cell_idx_t y,
                                                cell_idx_t z,
                                                Real_T & rho,
                                                Real_T & ux,
                                                Real_T & uy,
                                                Real_T & uz)
{
    const Real_T f0 = pdfField.get(x, y, z, 0);
    const Real_T f1 = pdfField.get(x, y, z, 1);
    const Real_T f2 = pdfField.get(x, y, z, 2);
    const Real_T f3 = pdfField.get(x, y, z, 3);
    const Real_T f4 = pdfField.get(x, y, z, 4);
    const Real_T f5 = pdfField.get(x, y, z, 5);
    const Real_T f6 = pdfField.get(x, y, z, 6);
    const Real_T f7 = pdfField.get(x, y, z, 7);
    const Real_T f8 = pdfField.get(x, y, z, 8);
    const Real_T f9 = pdfField.get(x, y, z, 9);
    const Real_T f10 = pdfField.get(x, y, z, 10);
    const Real_T f11 = pdfField.get(x, y, z, 11);
    const Real_T f12 = pdfField.get(x, y, z, 12);
    const Real_T f13 = pdfField.get(x, y, z, 13);
    const Real_T f14 = pdfField.get(x, y, z, 14);
    const Real_T f15 = pdfField.get(x, y, z, 15);
    const Real_T f16 = pdfField.get(x, y, z, 16);
    const Real_T f17 = pdfField.get(x, y, z, 17);
    const Real_T f18 = pdfField.get(x, y, z, 18);

    const Real_T vel0Term = f10 + f14 + f18 + f4 + f8;
    const Real_T momdensity0 = vel0Term - f13 - f17 - f3 - f7 - f9;
    const Real_T vel1Term = f11 + f15 + f7 + f1;
    const Real_T momdensity1 = vel1Term - f10 - f12 - f16 - f2 + f8 - f9;
    const Real_T vel2Term = f12 + f13 + f5;
    const Real_T deltaRho = vel0Term + vel1Term + vel2Term + f16 + f17 + f2 + f3 + f6 + f9 + f0;
    const Real_T momdensity2 = vel2Term + f11 + f14 - f15 - f16 - f17 - f18 - f6;

    const Real_T rawRho = deltaRho + Real_T(1);
    const Real_T safeRho = rawRho > Real_T(1e-12) ? rawRho : Real_T(1);
    const Real_T fx = forceField.get(x, y, z, 0);
    const Real_T fy = forceField.get(x, y, z, 1);
    const Real_T fz = forceField.get(x, y, z, 2);

    rho = safeRho;
    ux = momdensity0 / safeRho + Real_T(0.5) * fx / safeRho;
    uy = momdensity1 / safeRho + Real_T(0.5) * fy / safeRho;
    uz = momdensity2 / safeRho + Real_T(0.5) * fz / safeRho;
}

// Reason for edit start: replace the equilibrium-based top reconstruction with a free-slip mirror plus a relaxed density-sum correction so the custom top boundary behaves more like the generated slip wall while only nudging the incoming mass toward the previous top-cell state.
template <typename Real_T>
HOST_DEVICE_PREFIX void writeTopIncomingPopulationsWithSlipDensityCorrection(const GpuFieldView<Real_T> & pdfField,
                                                                             cell_idx_t x,
                                                                             cell_idx_t y,
                                                                             cell_idx_t zFluid,
                                                                             cell_idx_t zGhost)
{
    constexpr Real_T kDensityCorrectionRelaxation = Real_T(0.01); // Relaxation factor for nudging the incoming mass toward the previous top-cell state, can be tuned for stability vs responsiveness.

    const Real_T currentIncomingSum =
        pdfField.get(x, y, zFluid, 6) +
        pdfField.get(x, y, zFluid, 15) +
        pdfField.get(x, y, zFluid, 16) +
        pdfField.get(x, y, zFluid, 17) +
        pdfField.get(x, y, zFluid, 18);

    const Real_T mirroredIncomingSum =
        pdfField.get(x, y, zFluid, 5) +
        pdfField.get(x, y, zFluid, 11) +
        pdfField.get(x, y, zFluid, 12) +
        pdfField.get(x, y, zFluid, 13) +
        pdfField.get(x, y, zFluid, 14);

    const Real_T relaxedCorrection =
        kDensityCorrectionRelaxation * (currentIncomingSum - mirroredIncomingSum);
    const Real_T axisCorrection = relaxedCorrection * Real_T(1.0 / 3.0);
    const Real_T diagonalCorrection = relaxedCorrection * Real_T(1.0 / 6.0);

    pdfField.get(x, y, zGhost, 6) = pdfField.get(x, y, zFluid, 5) + axisCorrection;
    pdfField.get(x, y, zGhost, 15) = pdfField.get(x, y, zFluid, 11) + diagonalCorrection;
    pdfField.get(x, y, zGhost, 16) = pdfField.get(x, y, zFluid, 12) + diagonalCorrection;
    pdfField.get(x, y, zGhost, 17) = pdfField.get(x, y, zFluid, 13) + diagonalCorrection;
    pdfField.get(x, y, zGhost, 18) = pdfField.get(x, y, zFluid, 14) + diagonalCorrection;
}
// Reason for edit end: replace the equilibrium-based top reconstruction with a free-slip mirror plus a relaxed density-sum correction so the custom top boundary behaves more like the generated slip wall while only nudging the incoming mass toward the previous top-cell state.

template <typename Real_T>
GLOBAL_PREFIX void applyTopSlipZeroGradientKernel(GpuFieldView<Real_T> pdfField,
                                                  GpuFieldView<Real_T> forceField)
{
    const int32_t xGhostIndex = static_cast<int32_t>(blockIdx.x * blockDim.x + threadIdx.x);
    const int32_t yGhostIndex = static_cast<int32_t>(blockIdx.y * blockDim.y + threadIdx.y);
    const int32_t totalX = static_cast<int32_t>(pdfField.xSize + 2 * pdfField.ghostLayers);
    const int32_t totalY = static_cast<int32_t>(pdfField.ySize + 2 * pdfField.ghostLayers);

    if (xGhostIndex >= totalX || yGhostIndex >= totalY)
    {
        return;
    }

    const cell_idx_t x = static_cast<cell_idx_t>(xGhostIndex - pdfField.ghostLayers);
    const cell_idx_t y = static_cast<cell_idx_t>(yGhostIndex - pdfField.ghostLayers);
    const cell_idx_t zFluid = static_cast<cell_idx_t>(pdfField.zSize - 1);
    const cell_idx_t zGhost = static_cast<cell_idx_t>(pdfField.zSize);

    static_cast<void>(forceField);
    // Reason for edit start: switch the kernel from an equilibrium-based reconstruction to a reflected free-slip update with a gentle density-sum relaxation so the top boundary keeps slip-like momentum reflection without the strong turbulence inflation seen in the previous formulation.
    writeTopIncomingPopulationsWithSlipDensityCorrection(
        pdfField,
        x,
        y,
        zFluid,
        zGhost);
    // Reason for edit end: switch the kernel from an equilibrium-based reconstruction to a reflected free-slip update with a gentle density-sum relaxation so the top boundary keeps slip-like momentum reflection without the strong turbulence inflation seen in the previous formulation.
}

} // namespace internal

template <typename PdfField_T, typename ForceField_T>
class TopSlipZeroGradientBoundary
{
public:
    TopSlipZeroGradientBoundary(const std::shared_ptr<walberla::StructuredBlockForest> & blocks,
                                BlockDataID pdfFieldID,
                                BlockDataID forceFieldID,
                                uint_t domainNz,
                                bool enabled)
        : blocks_(blocks),
          pdfFieldID_(pdfFieldID),
          forceFieldID_(forceFieldID),
          domainNz_(domainNz),
          enabled_(enabled)
    {}

    bool isEnabled() const
    {
        return enabled_;
    }

    void operator()(walberla::IBlock * block, gpuStream_t stream = nullptr) const
    {
        run(block, stream);
    }

    void run(walberla::IBlock * block, gpuStream_t stream = nullptr) const
    {
        if (!enabled_)
        {
            return;
        }

        auto * pdfField = block->getData<PdfField_T>(pdfFieldID_);
        auto * forceField = block->getData<ForceField_T>(forceFieldID_);

        WALBERLA_ASSERT_NOT_NULLPTR(pdfField)
        WALBERLA_ASSERT_NOT_NULLPTR(forceField)

        walberla::Cell localOrigin{};
        localOrigin[0] = cell_idx_t(0);
        localOrigin[1] = cell_idx_t(0);
        localOrigin[2] = cell_idx_t(0);

        walberla::Cell globalOrigin{};
        blocks_->transformBlockLocalToGlobalCell(globalOrigin, *block, localOrigin);

        const cell_idx_t blockTopExclusive = static_cast<cell_idx_t>(globalOrigin[2] + pdfField->zSize());
        if (blockTopExclusive != static_cast<cell_idx_t>(domainNz_))
        {
            return;
        }

        const int totalX = static_cast<int>(pdfField->xSize() + 2 * pdfField->xOff());
        const int totalY = static_cast<int>(pdfField->ySize() + 2 * pdfField->yOff());
        const dim3 threads(internal::kTopSlipZeroGradientThreadsX, internal::kTopSlipZeroGradientThreadsY, 1u);
        const dim3 grid(
            static_cast<unsigned int>((totalX + static_cast<int>(threads.x) - 1) / static_cast<int>(threads.x)),
            static_cast<unsigned int>((totalY + static_cast<int>(threads.y) - 1) / static_cast<int>(threads.y)),
            1u);

        internal::applyTopSlipZeroGradientKernel<<<grid, threads, 0, stream>>>(
            internal::makeGpuFieldView(pdfField),
            internal::makeGpuFieldView(forceField));

        WALBERLA_GPU_CHECK(gpuPeekAtLastError());
    }

    std::function<void(walberla::IBlock *)> getSweep(gpuStream_t stream = nullptr) const
    {
        return [this, stream](walberla::IBlock * block) { this->run(block, stream); };
    }

private:
    std::shared_ptr<walberla::StructuredBlockForest> blocks_;
    BlockDataID pdfFieldID_{};
    BlockDataID forceFieldID_{};
    uint_t domainNz_{};
    bool enabled_{false};
};

} // namespace boundary
} // namespace turbine_core

#endif // TURBINECORE_ABLAPP_GPU_UNIFORM_TOPSLIPZEROGRADIENTBOUNDARY_H
// Reason for edit end: add a top-only sweep that keeps setup=Open and the bottom WFB path intact while approximating zero normal pressure gradient with slip-style kinematics at the top boundary.
