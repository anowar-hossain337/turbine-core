#pragma once

#ifndef TURBINECORE_ABLAPP_GPU_UNIFORM_TOPDAMPINGZONE_H
#define TURBINECORE_ABLAPP_GPU_UNIFORM_TOPDAMPINGZONE_H

#include <cstddef>
#include <functional>
#include <memory>

#include <blockforest/StructuredBlockForest.h>
#include <domain_decomposition/BlockDataID.h>
#include <domain_decomposition/IBlock.h>
#include <gpu/ErrorChecking.h>
#include <gpu/GPUField.h>
#include <gpu/GPUWrapper.h>

#include "walberla_helper/field/Field.h"

namespace turbine_core {
namespace damping {

namespace internal {

static constexpr int kThreadsX = 8;
static constexpr int kThreadsY = 8;
static constexpr int kThreadsZ = 4;

// This lightweight POD view is built on the host and passed into the kernel.
// It avoids constructing turbine_core::field::Field inside device-callable code,
// which was the source of the new host/device warnings from nvcc.
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

    HOST_DEVICE_PREFIX Real_T & get(cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f)
    {
        const std::ptrdiff_t idx =
            static_cast<std::ptrdiff_t>(f) * static_cast<std::ptrdiff_t>(fStride) +
            static_cast<std::ptrdiff_t>(x + ghostLayers) * static_cast<std::ptrdiff_t>(xStride) +
            static_cast<std::ptrdiff_t>(y + ghostLayers) * static_cast<std::ptrdiff_t>(yStride) +
            static_cast<std::ptrdiff_t>(z + ghostLayers) * static_cast<std::ptrdiff_t>(zStride);
        return data[idx];
    }
};

// Build the device-friendly view on the host once per launch.
template <typename Real_T>
HOST_PREFIX GpuFieldView<Real_T> makeGpuFieldView(walberla::gpu::GPUField<Real_T> * field)
{
    const cell_idx_t ghostLayers = field->xOff();

    return GpuFieldView<Real_T>{
        field->dataAt(-ghostLayers, -ghostLayers, -ghostLayers, 0),
        // Cast the waLBerla field metadata explicitly so nvcc does not report
        // aggregate-initializer narrowing while we populate the POD view.
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
GLOBAL_PREFIX void applyTopDampingKernel(GpuFieldView<Real_T> forceField,
                                         GpuFieldView<Real_T> velocityField,
                                         GpuFieldView<Real_T> meanVelocityField,
                                         int32_t globalZOffset,
                                         int32_t domainNz,
                                         int32_t dampingHeight,
                                         Real_T maxStrength,
                                         bool dampHorizontal,
                                         bool dampVertical)
{
    const int32_t x = static_cast<int32_t>(blockIdx.x * blockDim.x + threadIdx.x);
    const int32_t y = static_cast<int32_t>(blockIdx.y * blockDim.y + threadIdx.y);
    const int32_t z = static_cast<int32_t>(blockIdx.z * blockDim.z + threadIdx.z);

    // Access the POD members directly here; this view intentionally exposes raw
    // metadata instead of method calls so the kernel stays device-only.
    if (x >= static_cast<int32_t>(forceField.xSize) ||
        y >= static_cast<int32_t>(forceField.ySize) ||
        z >= static_cast<int32_t>(forceField.zSize))
    {
        return;
    }

    const int32_t cappedHeight = dampingHeight < domainNz ? dampingHeight : domainNz;
    const int32_t effectiveHeight = cappedHeight > 1 ? cappedHeight : 1;
    const int32_t dampingStart = domainNz - effectiveHeight;
    const int32_t globalZ = globalZOffset + z;

    if (globalZ < dampingStart)
    {
        return;
    }

    const Real_T eta = effectiveHeight <= 1
                           ? Real_T(1)
                           : Real_T(globalZ - dampingStart) / Real_T(effectiveHeight - 1);
    const Real_T sigma = maxStrength * eta * eta;

    if (dampHorizontal)
    {
        forceField.get(x, y, z, 0) += sigma * (meanVelocityField.get(x, y, z, 0) - velocityField.get(x, y, z, 0));
        forceField.get(x, y, z, 1) += sigma * (meanVelocityField.get(x, y, z, 1) - velocityField.get(x, y, z, 1));
    }

    if (dampVertical)
    {
        forceField.get(x, y, z, 2) += sigma * (-velocityField.get(x, y, z, 2));
    }
}

} // namespace internal


class TopDampingZone
{
public:
    TopDampingZone(const std::shared_ptr<walberla::StructuredBlockForest> & blocks,
                   BlockDataID forceFieldID,
                   BlockDataID velocityFieldID,
                   BlockDataID meanVelocityFieldID,
                   uint_t domainNz,
                   uint_t dampingHeight,
                   real_t maxStrength,
                   bool dampHorizontal,
                   bool dampVertical)
        : blocks_(blocks),
          forceFieldID_(forceFieldID),
          velocityFieldID_(velocityFieldID),
          meanVelocityFieldID_(meanVelocityFieldID),
          domainNz_(domainNz),
          dampingHeight_(dampingHeight),
          maxStrength_(maxStrength),
          dampHorizontal_(dampHorizontal),
          dampVertical_(dampVertical)
    {}

    bool isEnabled() const
    {
        return dampingHeight_ > uint_t(0) &&
               maxStrength_ > real_t(0) &&
               (dampHorizontal_ || dampVertical_);
    }

    void operator()(walberla::IBlock * block, gpuStream_t stream = nullptr) const
    {
        run(block, stream);
    }

    void run(walberla::IBlock * block, gpuStream_t stream = nullptr) const
    {
        if (!isEnabled())
        {
            return;
        }

        auto * forceField = block->getData<walberla::gpu::GPUField<real_t>>(forceFieldID_);
        auto * velocityField = block->getData<walberla::gpu::GPUField<real_t>>(velocityFieldID_);
        auto * meanVelocityField = block->getData<walberla::gpu::GPUField<real_t>>(meanVelocityFieldID_);

        WALBERLA_ASSERT_NOT_NULLPTR(forceField)
        WALBERLA_ASSERT_NOT_NULLPTR(velocityField)
        WALBERLA_ASSERT_NOT_NULLPTR(meanVelocityField)

        walberla::Cell localOrigin{};
        localOrigin[0] = cell_idx_t(0);
        localOrigin[1] = cell_idx_t(0);
        localOrigin[2] = cell_idx_t(0);

        walberla::Cell globalOrigin{};
        blocks_->transformBlockLocalToGlobalCell(globalOrigin, *block, localOrigin);

        const dim3 threads(internal::kThreadsX, internal::kThreadsY, internal::kThreadsZ);
        const dim3 grid(
            static_cast<unsigned int>((forceField->xSize() + threads.x - 1u) / threads.x),
            static_cast<unsigned int>((forceField->ySize() + threads.y - 1u) / threads.y),
            static_cast<unsigned int>((forceField->zSize() + threads.z - 1u) / threads.z));

        // Keep the old wrapper-based launch arguments as comments for traceability:
        // field::Field<real_t>(forceField),
        // field::Field<real_t>(velocityField),
        // field::Field<real_t>(meanVelocityField),
        // The new POD view removes the nvcc host/device warnings without changing
        // the actual damping formula.
        internal::applyTopDampingKernel<<<grid, threads, 0, stream>>>(
            internal::makeGpuFieldView(forceField),
            internal::makeGpuFieldView(velocityField),
            internal::makeGpuFieldView(meanVelocityField),
            static_cast<int32_t>(globalOrigin[2]),
            static_cast<int32_t>(domainNz_),
            static_cast<int32_t>(dampingHeight_),
            maxStrength_,
            dampHorizontal_,
            dampVertical_);

        WALBERLA_GPU_CHECK(gpuPeekAtLastError());
    }

    std::function<void(walberla::IBlock *)> getSweep(gpuStream_t stream = nullptr) const
    {
        return [this, stream](walberla::IBlock * block) { this->run(block, stream); };
    }

private:
    std::shared_ptr<walberla::StructuredBlockForest> blocks_;
    BlockDataID forceFieldID_;
    BlockDataID velocityFieldID_;
    BlockDataID meanVelocityFieldID_;
    uint_t domainNz_;
    uint_t dampingHeight_;
    real_t maxStrength_;
    bool dampHorizontal_;
    bool dampVertical_;
};

} // namespace damping
} // namespace turbine_core

#endif // TURBINECORE_ABLAPP_GPU_UNIFORM_TOPDAMPINGZONE_H
