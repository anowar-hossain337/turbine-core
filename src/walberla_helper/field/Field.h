
#pragma once
#ifndef TURBINECORE_FIELD_H
#define TURBINECORE_FIELD_H

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/math/Vector3.h"

#include <field/Layout.h>

namespace turbine_core {

    namespace field {

        //TODO assumes nrGhostLayers == 1
        template<typename Type_T>
        class Field {

        public:

            template< typename Field_T >
            HOST_DEVICE_PREFIX explicit Field( Field_T * field )
                    : size_{ field->xSize(), field->ySize(), field->zSize(), field->fSize() },
                      allocSize_{ field->xAllocSize(), field->yAllocSize(), field->zAllocSize(), field->fAllocSize() },
                      stride_{ field->xStride(), field->yStride(), field->zStride(), field->fStride() },
                      nGhostLayers_{ field->xOff() }
            {
                assert(field->layout() == walberla::field::fzyx && "Custom Field currently only supports fzyx layout!");
                data_ = field->dataAt(-nGhostLayers_, -nGhostLayers_, -nGhostLayers_, 0);
            }

            HOST_DEVICE_PREFIX auto xSize() const { return size_[0]; }
            HOST_DEVICE_PREFIX auto ySize() const { return size_[1]; }
            HOST_DEVICE_PREFIX auto zSize() const { return size_[2]; }
            HOST_DEVICE_PREFIX auto fSize() const { return size_[3]; }

            HOST_DEVICE_PREFIX auto nrOfGhostLayers() {
                return nGhostLayers_;
            }

            /**
             * \brief returns the value of the field for given local cell (x,y,z)
             *
             * If only global cell indices are given, use BlockInfo to obtain the local ones.
             */
            HOST_DEVICE_PREFIX       Type_T & get( cell_idx_t x, cell_idx_t y, cell_idx_t z ) {
                return this->get(x,y,z,0);
            }
            HOST_DEVICE_PREFIX const Type_T & get( cell_idx_t x, cell_idx_t y, cell_idx_t z ) const {
                return this->get(x,y,z,0);
            }
            HOST_DEVICE_PREFIX       Type_T & get( cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f ) {
                const Field<Type_T> & const_this = *this;
                return const_cast<Type_T&>(const_this.get(x,y,z,f));
            }
            HOST_DEVICE_PREFIX const Type_T & get( cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f ) const {
                //TODO validate coordinates

                const uint_t idx = f * stride_[3] + (x+nGhostLayers_) * stride_[0] + (y+nGhostLayers_) * stride_[1] + (z+nGhostLayers_) * stride_[2];

                //TODO check if index is valid

                return *(data_ + idx);
            }

            HOST_DEVICE_PREFIX const Type_T & get( const Vector3<cell_idx_t> & xyz ) const {
                return get(xyz[0], xyz[1], xyz[2], 0);
            }

            HOST_DEVICE_PREFIX       Type_T & get( const Vector3<cell_idx_t> & xyz  ) {
                return get(xyz[0], xyz[1], xyz[2], 0);
            }

            HOST_DEVICE_PREFIX const Type_T & get( const Vector3<cell_idx_t> & xyz, const cell_idx_t f ) const {
                return get(xyz[0], xyz[1], xyz[2], f);
            }

            HOST_DEVICE_PREFIX       Type_T & get( const Vector3<cell_idx_t> & xyz, const cell_idx_t f ) {
                return get(xyz[0], xyz[1], xyz[2], f);
            }

            HOST_DEVICE_PREFIX const Type_T & get( const Vector3<cell_idx_t> & xyz, const uint_t f ) const {
                return get(xyz[0], xyz[1], xyz[2], f);
            }

            HOST_DEVICE_PREFIX       Type_T & get( const Vector3<cell_idx_t> & xyz, const uint_t f ) {
                return get(xyz[0], xyz[1], xyz[2], f);
            }

            HOST_DEVICE_PREFIX const Type_T & getF( Type_T * xyz0, uint_t f ) const {
                assert(f < size_[3] && "Cannot access f that exceeds fSize of field.");
                return *(xyz0 + f * stride_[3]);
            }

            HOST_DEVICE_PREFIX Type_T & getF( Type_T * xyz0, uint_t f ) {
                assert(f < size_[3] && "Cannot access f that exceeds fSize of field.");
                return *(xyz0 + f * stride_[3]);
            }

            HOST_DEVICE_PREFIX const Type_T & getF( Type_T * xyz0, cell_idx_t f ) const {
                assert(f < size_[3] && "Cannot access f that exceeds fSize of field.");
                return *(xyz0 + f * stride_[3]);
            }

            HOST_DEVICE_PREFIX Type_T & getF( Type_T * xyz0, cell_idx_t f ) {
                assert(f < size_[3] && "Cannot access f that exceeds fSize of field.");
                return *(xyz0 + f * stride_[3]);
            }

        private:

            // stored as x,y,z,f
            const uint_t size_[4];
            const uint_t allocSize_[4];
            const uint_t stride_[4];

            const uint_t nGhostLayers_;

            Type_T * RESTRICT data_;

        };

    }

}

#endif //TURBINECORE_FIELD_H
