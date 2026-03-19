
#pragma once

#ifndef TURBINECORE_TRILINEARFIELDINTERPOLATOR_H
#define TURBINECORE_TRILINEARFIELDINTERPOLATOR_H

#include <cassert>

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

#include "wind_turbine_core/math/Vector3.h"

#include "walberla_helper/field/MacroscopicVariableCalculator.h"

#include "walberla_helper/blockforest/BlockInfo.h"

namespace turbine_core {

    namespace projectors {

        template< typename Field_T, typename InterpolatorType_T >
        class TrilinearFieldInterpolatorBase {

        public:

            HOST_DEVICE_PREFIX TrilinearFieldInterpolatorBase( const blockforest::BlockInfo & blockInfo,
                                                               Field_T * baseField )
                    : blockInfo_(&blockInfo), baseField_(baseField)
            {
                assert(baseField->nrOfGhostLayers() > uint_t(0) && "field for tri-linear field interpolator needs at least one ghost layer");
            }

            template<typename ForwardIterator_T>
            HOST_DEVICE_PREFIX void get( const Vector3<real_t> & position, ForwardIterator_T interpolationResultBegin ) {
                get( position[0], position[1], position[2], interpolationResultBegin );
            }

            template< typename ForwardIterator_T >
            HOST_DEVICE_PREFIX void get( const real_t x, const real_t y, const real_t z, ForwardIterator_T interpolationResultBegin ) {

                auto blockBB = blockInfo_->getAABB();

                assert(blockBB.contains(x,y,z) && "Interpolation position must be contained in block AABB.");

                auto containingCell = blockInfo_->getBlockLocalCell(x,y,z);
                auto containingCellCenter = blockInfo_->getBlockLocalCellCenter(containingCell);

                const auto dx = blockInfo_->dx();
                const auto dy = blockInfo_->dy();
                const auto dz = blockInfo_->dz();

                const auto xNeighbor = cell_idx_t( floor( (x - containingCellCenter[0]) / dx ) );
                const auto yNeighbor = cell_idx_t( floor( (y - containingCellCenter[1]) / dy ) );
                const auto zNeighbor = cell_idx_t( floor( (z - containingCellCenter[2]) / dz ) );

                // define the 8 nearest cells required for the tri-linear interpolation
                // the cell 'ccc' is the one with the smallest x-, y-, and z-indices

                Vector3<cell_idx_t> cells[8];
                for(uint_t k = 0; k < 2; ++k)
                    for(uint_t j = 0; j < 2; ++j)
                        for(uint_t i = 0; i < 2; ++i)
                            cells[4*k+2*j+i] = Vector3<cell_idx_t>(
                                    containingCell[0] + xNeighbor + cell_idx_t(i),
                                    containingCell[1] + yNeighbor + cell_idx_t(j),
                                    containingCell[2] + zNeighbor + cell_idx_t(k)
                            );

                // weighting
                const real_t inv_totalVolume = real_t(1) / ( blockInfo_->dx() * blockInfo_->dy() * blockInfo_->dz() );
                auto cccCellCenter = blockInfo_->getBlockLocalCellCenter(cells[0]);

                const real_t dxc[] = {cccCellCenter[0] + dx - x, x - cccCellCenter[0]};
                const real_t dyc[] = {cccCellCenter[1] + dy - y, y - cccCellCenter[1]};
                const real_t dzc[] = {cccCellCenter[2] + dz - z, z - cccCellCenter[2]};

                real_t weighting[8];
                for(uint_t k = 0; k < 2; ++k)
                    for(uint_t j = 0; j < 2; ++j)
                        for(uint_t i = 0; i < 2; ++i)
                            weighting[4*k+2*j+i] = dxc[i] * dyc[j] * dzc[k] * inv_totalVolume;

                auto derived = static_cast<InterpolatorType_T&>(*this);

                for(uint_t i = 0; i < 8; ++i) {
		        //TODO does not account for ghost layers in non-BC regions...
//                    if(!blockBB.contains(cells[i]))
//                        cells[i] = containingCell;
                    derived.addWeightedCellValue( interpolationResultBegin, cells[i], weighting[i] );
                }

            }

        protected:

            const blockforest::BlockInfo * blockInfo_;

            Field_T * baseField_;

        };


        template<typename Field_T>
        class TrilinearFieldInterpolator : public TrilinearFieldInterpolatorBase<Field_T, TrilinearFieldInterpolator<Field_T>> {

            friend TrilinearFieldInterpolatorBase<Field_T, TrilinearFieldInterpolator<Field_T>>;

        public:

            HOST_DEVICE_PREFIX TrilinearFieldInterpolator( const blockforest::BlockInfo & blockInfo,
                                                           Field_T * baseField )
                    : TrilinearFieldInterpolatorBase<Field_T,TrilinearFieldInterpolator<Field_T>>(blockInfo, baseField)
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

        template<typename Stencil_T, typename PdfField_T, typename ForceField_T, bool zeroCentering, bool compressible>
        class TrilinearMacroscopicFieldInterpolator : public TrilinearFieldInterpolatorBase<PdfField_T, TrilinearMacroscopicFieldInterpolator<Stencil_T, PdfField_T, ForceField_T, zeroCentering, compressible>> {

            using OwnType_T = TrilinearMacroscopicFieldInterpolator<Stencil_T, PdfField_T, ForceField_T, zeroCentering, compressible>;

            friend TrilinearFieldInterpolatorBase<PdfField_T, OwnType_T>;

        public:

            HOST_DEVICE_PREFIX TrilinearMacroscopicFieldInterpolator( const blockforest::BlockInfo & blockInfo,
                                                                      PdfField_T * pdfField, ForceField_T * forceField )
                    : TrilinearFieldInterpolatorBase<PdfField_T,OwnType_T>(blockInfo, pdfField)
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

#endif //TURBINECORE_TRILINEARFIELDINTERPOLATOR_H
