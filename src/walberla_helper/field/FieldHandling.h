
#pragma once
#ifndef TURBINECORE_FIELDHANDLING_H
#define TURBINECORE_FIELDHANDLING_H

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/math/Vector3.h"

#include <field/Layout.h>
#include <field/allocation/FieldAllocator.h>

namespace turbine_core {

    namespace field {

        // data handling for loading a field of type FlattenedField_T from file
        template< typename FlattenedField_T >
        class FlattenedFieldHandling : public walberla::field::BlockDataHandling< FlattenedField_T >
        {
        public:

            using Value_T = typename FlattenedField_T::value_type;

            FlattenedFieldHandling(const std::weak_ptr< walberla::StructuredBlockStorage >& blocks, const uint_t numberGhostLayer,
                                   const walberla::field::Layout & layout,
                                   const std::shared_ptr<walberla::field::FieldAllocator<Value_T> > &alloc = std::shared_ptr<walberla::field::FieldAllocator<Value_T>>())
                    : blocks_(blocks), alloc_(alloc), numberGhostLayer_(numberGhostLayer), layout_(layout)
            {}

        protected:
            FlattenedField_T* allocate(walberla::IBlock* const block) override { return allocateDispatch(block); }

            FlattenedField_T* reallocate(walberla::IBlock* const block) override { return allocateDispatch(block); }

        private:
            std::weak_ptr< walberla::StructuredBlockStorage > blocks_;

            const uint_t numberGhostLayer_;
            const walberla::field::Layout layout_;
            const std::shared_ptr<walberla::field::FieldAllocator<Value_T> > alloc_;

            FlattenedField_T* allocateDispatch(walberla::IBlock* const block)
            {
                WALBERLA_ASSERT_NOT_NULLPTR(block)

                auto blocks = blocks_.lock();
                WALBERLA_CHECK_NOT_NULLPTR(blocks)

                return new FlattenedField_T(blocks->getNumberOfXCells(*block), blocks->getNumberOfYCells(*block),
                                            blocks->getNumberOfZCells(*block), numberGhostLayer_, layout_, alloc_);
            }
        }; // class FlattenedFieldHandling

    }

}

#endif //TURBINECORE_FIELDHANDLING_H
