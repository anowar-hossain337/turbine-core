#pragma once

#ifndef TURBINECORE_WALBERLAHELPER_FIELD_ALL_H
#define TURBINECORE_WALBERLAHELPER_FIELD_ALL_H

#include "walberla_helper/field/distributors/KernelFieldDistributor.h"
#include "walberla_helper/field/distributors/NearestNeighbourFieldDistributor.h"

#include "walberla_helper/field/interpolators/KernelFieldInterpolator.h"
#include "walberla_helper/field/interpolators/NearestNeighbourFieldInterpolator.h"
#include "walberla_helper/field/interpolators/TrilinearFieldInterpolator.h"

#include "walberla_helper/field/projector_kernels/GaussianKernel.h"
#include "walberla_helper/field/projector_kernels/M4ProjectionKernel.h"
#include "walberla_helper/field/projector_kernels/SmoothedDiracDeltaKernel.h"

#include "walberla_helper/field/FieldCopy.h"
#include "walberla_helper/field/FieldHandling.h"

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
    #include "walberla_helper/field/FieldCopyGPU.h"
#endif


#endif //TURBINECORE_WALBERLAHELPER_FIELD_ALL_H
