
#pragma once

#ifndef TURBINECORE_BLOCKINFO_H
#define TURBINECORE_BLOCKINFO_H

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

#include "wind_turbine_core/math/Vector3.h"
#include "wind_turbine_core/math/AABB.h"

#include <memory>

#include <domain_decomposition/IBlock.h>
#include <blockforest/StructuredBlockForest.h>

namespace turbine_core {

    namespace blockforest {

        struct BlockInfo {

        public:

            HOST_PREFIX BlockInfo( walberla::IBlock * const block, walberla::blockforest::StructuredBlockForest * const forest )
                    : level_(forest->getLevel(*block)),
                      dx_{forest->dx(level_),
                          forest->dy(level_),
                          forest->dz(level_)},
                      blockBB_(block->getAABB())
            {}

            HOST_PREFIX BlockInfo( walberla::IBlock * const block, const std::shared_ptr<walberla::blockforest::StructuredBlockForest> & forest )
                    : BlockInfo(block,forest.get())
            {}

            HOST_DEVICE_PREFIX BlockInfo( const BlockInfo & blockInfo )
                    : level_(blockInfo.level_),
                      dx_{blockInfo.dx_[0], blockInfo.dx_[1], blockInfo.dx_[2]},
                      blockBB_(blockInfo.blockBB_)

            {}

            HOST_DEVICE_PREFIX auto dx() const { return dx_[0]; }
            HOST_DEVICE_PREFIX auto dy() const { return dx_[1]; }
            HOST_DEVICE_PREFIX auto dz() const { return dx_[2]; }

            HOST_DEVICE_PREFIX auto & getAABB() const {
                return blockBB_;
            }

            HOST_DEVICE_PREFIX Vector3<cell_idx_t> getBlockLocalCell( const real_t x, const real_t y, const real_t z ) const {

                Vector3<cell_idx_t> localCell;

                localCell[0] = cell_idx_t( floor( ( x - blockBB_.xMin() ) / dx_[0] ) );
                localCell[1] = cell_idx_t( floor( ( y - blockBB_.yMin() ) / dx_[1] ) );
                localCell[2] = cell_idx_t( floor( ( z - blockBB_.zMin() ) / dx_[2] ) );

                return localCell;
            }

            HOST_DEVICE_PREFIX Vector3<cell_idx_t> getBlockLocalCell( const Vector3<real_t> & position ) const {
                return getBlockLocalCell(position[0], position[1], position[2]);
            }


            HOST_DEVICE_PREFIX Vector3<real_t> getBlockLocalCellCenter( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const {

                real_t centerX = blockBB_.xMin() + ( real_t(x) + real_t(0.5) ) * dx_[0];
                real_t centerY = blockBB_.yMin() + ( real_t(y) + real_t(0.5) ) * dx_[1];
                real_t centerZ = blockBB_.zMin() + ( real_t(z) + real_t(0.5) ) * dx_[2];

                return {centerX,centerY,centerZ};
            }

            HOST_DEVICE_PREFIX Vector3<real_t> getBlockLocalCellCenter( const Vector3<cell_idx_t> & localCellIndices ) const {
                return getBlockLocalCellCenter(localCellIndices[0], localCellIndices[1], localCellIndices[2]);
            }


        private:

            const uint_t level_;
            const real_t dx_[3];

            const AABB blockBB_;

        };

    }

}

#endif //TURBINECORE_BLOCKINFO_H
