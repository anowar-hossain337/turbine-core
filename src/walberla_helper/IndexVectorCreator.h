
#pragma once

#ifndef TURBINECORE_INDEXVECTORCREATOR_H
#define TURBINECORE_INDEXVECTORCREATOR_H

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/math/Vector3.h"

#include <vector>

#include <blockforest/StructuredBlockForest.h>
#include <field/FlagField.h>

namespace turbine_core {

    namespace field {

        struct IndexInfo {
            cell_idx_t x;
            cell_idx_t y;
            cell_idx_t z;

            HOST_DEVICE_PREFIX IndexInfo( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
                    : x(x), y(y), z(z)
            {}

            HOST_DEVICE_PREFIX bool operator==( const IndexInfo & info ) const {
                return (x == info.x) && (y == info.y) && (z == info.z);
            }

            template<typename CellIndex_T>
            HOST_DEVICE_PREFIX bool operator==( const Vector3<CellIndex_T> & indices ) {
                return (x == static_cast<cell_idx_t>(indices[0])) &&
                       (y == static_cast<cell_idx_t>(indices[1])) &&
                       (z == static_cast<cell_idx_t>(indices[2]));
            }
        };

        class IndexVector {

        public:

            HOST_PREFIX ~IndexVector() {
#ifdef WALBERLA_BUILD_WITH_CUDA
                cudaFree( gpuVector_ );
#endif
            }

            HOST_PREFIX void syncGPU() {
#ifdef WALBERLA_BUILD_WITH_CUDA
                cudaFree(gpuVector_);
                cudaMalloc((void**)&gpuVector_, cpuVector_.size() * sizeof(IndexInfo));
                cudaMemcpy( gpuVector_, cpuVector_.data(), cpuVector_.size() * sizeof(IndexInfo), cudaMemcpyHostToDevice );
#endif
            }

            HOST_PREFIX auto & cpuInfoVector() {
                return cpuVector_;
            }

            HOST_PREFIX auto * cpuPointer() {
                return cpuVector_.data();
            }

            HOST_PREFIX auto * gpuPointer() {
                return gpuVector_;
            }

            HOST_PREFIX uint_t nIndexInfo() const {
                return cpuVector_.size();
            }

            HOST_PREFIX bool operator==( const IndexVector & idx ) {
                return cpuVector_ == idx.cpuVector_;
            }

        private:

            std::vector<IndexInfo> cpuVector_{};
            IndexInfo * gpuVector_{};

        };

        class IndexVectorCreator {

        public:

            HOST_PREFIX explicit IndexVectorCreator( const std::shared_ptr<walberla::blockforest::StructuredBlockForest> & forest )
            : forest_(forest)
            {
                auto createIndexVectors = []( walberla::IBlock * const, walberla::domain_decomposition::StructuredBlockStorage * const ) {return new IndexVector();};
                indexVectorID_ = forest_->addStructuredBlockData<IndexVector>( createIndexVectors, "IndexField" );
            }

            HOST_PREFIX auto indexVectorID() const {
                return indexVectorID_;
            }

            template< typename FlagField_T >
            HOST_PREFIX void fillFromFlagField( const walberla::BlockDataID & flagFieldID, const walberla::field::FlagUID & domainUID) {

                for( auto blockIt = forest_->begin(); blockIt != forest_->end(); ++blockIt ) {

                    auto * indexVector = blockIt->getData<IndexVector>( indexVectorID_ );
                    auto & cpuIndexVector = indexVector->cpuInfoVector();

                    auto * flagField   = blockIt->getData<FlagField_T>( flagFieldID );

                    if( ! flagField->flagExists(domainUID) )
                        return;

                    auto domainFlag = flagField->getFlag(domainUID);

                    cpuIndexVector.clear();

                    // iterate over field
                    for( auto it = flagField->begin(); it != flagField->end(); ++it ) {

                        // store non-domain cells
                        if( ! walberla::isFlagSet( it, domainFlag ) ) {
                            IndexInfo info {it.x(), it.y(), it.z()};
                            cpuIndexVector.push_back(info);
                        }

                    }

                    // copy cpu data to gpu
                    indexVector->syncGPU();

                }

            }

        private:

            walberla::BlockDataID indexVectorID_;

            const std::shared_ptr<walberla::StructuredBlockForest> forest_{};

        };

    }

}

#endif //TURBINECORE_INDEXVECTORCREATOR_H
