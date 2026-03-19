
#pragma once

#ifndef TURBINECORE_NEARESTNEIGHBOURFIELDINTERPOLATOR_H
#define TURBINECORE_NEARESTNEIGHBOURFIELDINTERPOLATOR_H

#include <cassert>

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

#include "wind_turbine_core/math/Vector3.h"

#include "walberla_helper/blockforest/BlockInfo.h"

namespace turbine_core {

    namespace projectors {

        template< typename Field_T >
        class NearestNeighbourFieldInterpolator {

        public:

            HOST_DEVICE_PREFIX NearestNeighbourFieldInterpolator( const blockforest::BlockInfo & blockInfo,
                                                                  Field_T * baseField )
                    : blockInfo_(&blockInfo), baseField_(baseField)
            {}

            template<typename ForwardIterator_T>
            HOST_DEVICE_PREFIX void get( const Vector3<real_t> & position, ForwardIterator_T interpolationResultBegin ) {
                get( position[0], position[1], position[2], interpolationResultBegin );
            }

            template< typename ForwardIterator_T >
            HOST_DEVICE_PREFIX void get( const real_t & x, const real_t & y, const real_t & z, ForwardIterator_T interpolationResultBegin ) {

                assert(blockInfo_->getAABB().contains(x,y,z) && "Interpolation position must be contains in block AABB.");

                auto containingCell = blockInfo_->getBlockLocalCell(x,y,z);
                auto containingCellCenter = blockInfo_->getBlockLocalCellCenter(containingCell);

                bool isBoundaryCell = false;
//                for( uint_t i = 0; i < nIndexInfo_; ++i ) {
//                    if(indexInfo_[i] == containingCell) {
//                        isBoundaryCell = true;
//                        break;
//                    }
//                }

                if (!isBoundaryCell) {
                    addCellValue(interpolationResultBegin, containingCell);
                    return;
                } else {
                    const auto xNeighbor = cell_idx_t( floor( x - containingCellCenter[0] ) );
                    const auto yNeighbor = cell_idx_t( floor( y - containingCellCenter[1] ) );
                    const auto zNeighbor = cell_idx_t( floor( z - containingCellCenter[2] ) );

                    const cell_idx_t xMin = containingCell[0] + xNeighbor;
                    const cell_idx_t yMin = containingCell[1] + yNeighbor;
                    const cell_idx_t zMin = containingCell[2] + zNeighbor;

                    for( cell_idx_t zC = zMin; zC <= zMin + cell_idx_t(1); ++zC) {
                        for( cell_idx_t yC = yMin; yC <= yMin + cell_idx_t(1); ++yC) {
                            for( cell_idx_t xC = xMin; xC <= xMin + cell_idx_t(1); ++xC) {
                                Vector3<cell_idx_t> curCell(xC,yC,zC);
                                if( blockInfo_->getAABB().contains(real_t(xC), real_t(yC), real_t(zC)) ) {

                                    bool isDomainCell = true;
//                                    for(uint_t i = 0; i < nIndexInfo_; ++i) {
//                                        if(indexInfo_[i] == curCell) {
//                                            isDomainCell = false;
//                                            break;
//                                        }
//                                    }

                                    if (isDomainCell) {
                                        addCellValue(interpolationResultBegin, curCell);
                                        return;
                                    }
                                }
                            }
                        }
                    }
                }
            }

        private:

            template< typename ForwardIterator_T >
            HOST_DEVICE_PREFIX void addCellValue(ForwardIterator_T interpolationResultBegin, const Vector3<cell_idx_t> & curCell)
            {
                for (uint_t f = 0; f < baseField_->fSize(); ++f) {
                    auto value = baseField_->get(curCell, f);
                    assert(value == value && "NaN found in component when interpolating from cell.\n");
                    *interpolationResultBegin = value;
                    ++interpolationResultBegin;
                }
            }

            const blockforest::BlockInfo * blockInfo_;

            Field_T * baseField_;
        };

    }
}

#endif //TURBINECORE_NEARESTNEIGHBOURFIELDINTERPOLATOR_H
