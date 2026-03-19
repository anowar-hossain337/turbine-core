
#pragma once

#ifndef TURBINECORE_KERNELFIELDINTERPOLATOR_H
#define TURBINECORE_KERNELFIELDINTERPOLATOR_H

#include <cassert>

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

#include "wind_turbine_core/math/Vector3.h"

#include "walberla_helper/field/MacroscopicVariableCalculator.h"

#include "walberla_helper/blockforest/BlockInfo.h"

namespace turbine_core {

    namespace projectors {

        template< typename Field_T, typename InterpolatorType_T, typename Kernel_T >
        class KernelFieldInterpolatorBase {

        public:

            using Kernel = Kernel_T;

            HOST_DEVICE_PREFIX KernelFieldInterpolatorBase( const blockforest::BlockInfo & blockInfo,
                                                            Field_T * baseField )
                    : blockInfo_(&blockInfo), baseField_(baseField)
            {
                assert( (baseField->nrOfGhostLayers() >= max(max(static_cast<uint_t>(Kernel_T::xWidth), static_cast<uint_t>(Kernel_T::yWidth)), static_cast<uint_t>(Kernel_T::zWidth)))
                       && "field for kernel interpolator needs sufficiently many ghost layer");
            }

            template<typename ForwardIterator_T>
            HOST_DEVICE_PREFIX void get( const Vector3<real_t> & position, ForwardIterator_T interpolationResultBegin ) {
                get( position[0], position[1], position[2], interpolationResultBegin );
            }

            template< typename ForwardIterator_T >
            HOST_DEVICE_PREFIX void get( const real_t x, const real_t y, const real_t z, ForwardIterator_T interpolationResultBegin ) {

                auto centerCell = blockInfo_->getBlockLocalCell( x, y, z );

                const real_t dx = blockInfo_->dx();
                const real_t dy = blockInfo_->dy();
                const real_t dz = blockInfo_->dz();

                const Vector3<uint_t> neighborhoodSize {uint_t(Kernel_T::xWidth), uint_t(Kernel_T::yWidth), uint_t(Kernel_T::zWidth)};

                math::GenericAABB<cell_idx_t> cellNeighborhood( centerCell[0] - cell_idx_t(neighborhoodSize[0]), centerCell[1] - cell_idx_t(neighborhoodSize[1]), centerCell[2] - cell_idx_t(neighborhoodSize[2]),
                                                                centerCell[0] + cell_idx_t(neighborhoodSize[0]), centerCell[1] + cell_idx_t(neighborhoodSize[1]), centerCell[2] + cell_idx_t(neighborhoodSize[2]) );


                uint_t counter{0};
                real_t sumOfWeights {0};

                // get distribution weights and count available cells in surrounding cells
                for( cell_idx_t zC = cellNeighborhood.zMin(); zC <= cellNeighborhood.zMax(); ++zC ) {
                    for (cell_idx_t yC = cellNeighborhood.yMin(); yC <= cellNeighborhood.yMax(); ++yC) {
                        for (cell_idx_t xC = cellNeighborhood.xMin(); xC <= cellNeighborhood.xMax(); ++xC) {

                            Vector3<cell_idx_t> curCell{xC,yC,zC};
                            Vector3<real_t> curCellCenter = blockInfo_->getBlockLocalCellCenter(curCell);

                            weights_[counter] = Kernel_T::getWeight(x, y, z, curCellCenter[0], curCellCenter[1], curCellCenter[2],
                                                                    dx, dy, dz);
                            sumOfWeights += weights_[counter];
                            ++counter;
                        }
                    }
                }

                const real_t scalingFactor = real_t(1) / sumOfWeights;

                auto derived = static_cast<InterpolatorType_T&>(*this);

                counter = uint_t(0);
                for( cell_idx_t zC = cellNeighborhood.zMin(); zC <= cellNeighborhood.zMax(); ++zC ) {
                    for (cell_idx_t yC = cellNeighborhood.yMin(); yC <= cellNeighborhood.yMax(); ++yC) {
                        for (cell_idx_t xC = cellNeighborhood.xMin(); xC <= cellNeighborhood.xMax(); ++xC) {
                            if (weights_[counter] > real_t(0)) {
                                derived.addWeightedCellValue( interpolationResultBegin, Vector3<cell_idx_t>(xC,yC,zC), scalingFactor * weights_[counter] );
                            }
                            ++counter;
                        }
                    }
                }

            }

        protected:

            const blockforest::BlockInfo * blockInfo_;

            Field_T * baseField_;

            real_t weights_[(uint_t(1) + uint_t(2) * static_cast<uint_t>(Kernel_T::xWidth)) *
                            (uint_t(1) + uint_t(2) * static_cast<uint_t>(Kernel_T::yWidth)) *
                            (uint_t(1) + uint_t(2) * static_cast<uint_t>(Kernel_T::zWidth)) ]{};

        };


        template<typename Field_T, typename Kernel_T>
        class KernelFieldInterpolator : public KernelFieldInterpolatorBase<Field_T, KernelFieldInterpolator<Field_T, Kernel_T>, Kernel_T> {

            friend KernelFieldInterpolatorBase<Field_T, KernelFieldInterpolator<Field_T, Kernel_T>, Kernel_T>;

        public:

            HOST_DEVICE_PREFIX KernelFieldInterpolator( const blockforest::BlockInfo & blockInfo,
                                                        Field_T * baseField )
                    : KernelFieldInterpolatorBase<Field_T,KernelFieldInterpolator<Field_T, Kernel_T>, Kernel_T>(blockInfo, baseField)
            {}

        private:

            template< typename ForwardIterator_T >
            HOST_DEVICE_PREFIX void addWeightedCellValue( ForwardIterator_T interpolationResultBegin, const Vector3<cell_idx_t> & curCell, const real_t & weighting )
            {
                for( uint_t f = 0; f < this->baseField_->fSize(); ++f ) {

                    auto value = this->baseField_->get( curCell, f);

                    assert( (value == value) && "NaN found in component when interpolating from cell." );

                    *interpolationResultBegin += weighting * value;
                    ++interpolationResultBegin;
                }
            }

        };

        template<typename Kernel_T, typename Stencil_T, typename PdfField_T, typename ForceField_T, bool zeroCentering, bool compressible>
        class KernelMacroscopicFieldInterpolator : public KernelFieldInterpolatorBase<PdfField_T, KernelMacroscopicFieldInterpolator<Kernel_T, Stencil_T, PdfField_T, ForceField_T, zeroCentering, compressible>, Kernel_T> {

            using OwnType_T = KernelMacroscopicFieldInterpolator<Kernel_T, Stencil_T, PdfField_T, ForceField_T, zeroCentering, compressible>;

            friend KernelFieldInterpolatorBase<PdfField_T, OwnType_T, Kernel_T>;

        public:

            HOST_DEVICE_PREFIX KernelMacroscopicFieldInterpolator( const blockforest::BlockInfo & blockInfo,
                                                                   PdfField_T * pdfField, ForceField_T * forceField )
                    : KernelFieldInterpolatorBase<PdfField_T,OwnType_T, Kernel_T>(blockInfo, pdfField)
            {
                calculator_.setPdfField(pdfField);
                calculator_.setForceField(forceField);
            }

        private:

            template< typename ForwardIterator_T >
            HOST_DEVICE_PREFIX void addWeightedCellValue( ForwardIterator_T interpolationResultBegin, const Vector3<cell_idx_t> & curCell, const real_t & weighting )
            {
                real_t density{};
                Vector3<real_t> velocity{};
                calculator_.get(curCell, density, velocity);

                assert( (density == density) && "NaN found in density when interpolating from cell." );
                assert( (velocity[0] == velocity[0]) && "NaN found in velocity[0] when interpolating from cell." );
                assert( (velocity[1] == velocity[1]) && "NaN found in velocity[1] when interpolating from cell." );
                assert( (velocity[2] == velocity[2]) && "NaN found in velocity[2] when interpolating from cell." );

                (*interpolationResultBegin) += weighting * density;
                ++interpolationResultBegin;

                for( uint_t f = 0; f < 3; ++f ) {
                    *interpolationResultBegin += weighting * velocity[f];
                    ++interpolationResultBegin;
                }
            }

            field::DensityAndVelocityCalculator<Stencil_T,PdfField_T,ForceField_T,zeroCentering,compressible> calculator_;

        };

    }
}

#endif //TURBINECORE_KERNELFIELDINTERPOLATOR_H
