
#pragma once

#ifndef TURBINECORE_CUSTOMMEMORYBUFFER_H
#define TURBINECORE_CUSTOMMEMORYBUFFER_H

namespace turbine_core {

    namespace mpi {


        //*******************************************************************************************************************
        /*!
         * Simple buffer class that supports memory allocators, e.g. for pinned host memory or GPU memory
         *
         * \ingroup cuda
         *
         * In contrast to core::mpi::Buffer this class does not support stream operators "<<" and ">>" because these
         * operators imply serial (un)packing which is not feasible on the GPU.
         * The allocator template has to provide:
         *   - static void *allocate( size_t size )
         *   - void deallocate( void *ptr )
         *   - void memcpy( void *dst, void *src, size_t count )
         *
         * The buffer has a beginning, a current position and an end position. Here is an overview of the most important
         * operations:
         *   - clear: reset current position to begin, does not change size
         *   - advance: moves current position number of bytes forward and returns pointer to the old current position
         *              two versions are available, one that automatically resizes and reallocates the buffer, and one that
         *              fails if not enough space is available
         */
        //*******************************************************************************************************************
        template<typename Allocator>
        class CustomMemoryBuffer
        {
        public:
            typedef uint8_t ElementType;

            explicit CustomMemoryBuffer();
            explicit CustomMemoryBuffer( std::size_t initSize );
            CustomMemoryBuffer( const CustomMemoryBuffer &pb );
            ~CustomMemoryBuffer();
            CustomMemoryBuffer &operator=( const CustomMemoryBuffer &pb );

            void resize( std::size_t newSize );
            inline std::size_t allocSize() const { return std::size_t(end_ - begin_); }
            inline std::size_t size() const { return std::size_t(cur_ - begin_); }
            ElementType *ptr() const { return begin_; }
            ElementType *cur() const { return cur_; }

            inline bool isEmpty() const {return (size() == std::size_t(0));}

            inline void clear() { cur_ = begin_; }

            ElementType *advance( std::size_t bytes );
            ElementType *advanceNoResize( std::size_t bytes );

            template<typename T>
            T *advance( std::size_t bytes ) { return reinterpret_cast<T *>( advance( bytes * sizeof( T ))); }
            template<typename T>
            T *advanceNoResize( std::size_t bytes ) { return reinterpret_cast<T *>( advanceNoResize( bytes * sizeof( T ))); }

        private:
            ElementType *begin_ = nullptr;
            ElementType *cur_ = nullptr;
            ElementType *end_ = nullptr;
        };

        template<typename Allocator>
        CustomMemoryBuffer<Allocator>::CustomMemoryBuffer() = default;

        template<typename Allocator>
        CustomMemoryBuffer<Allocator>::CustomMemoryBuffer( std::size_t initSize )
        {
            if( initSize > 0 )
            {
                begin_ = reinterpret_cast<ElementType *>(Allocator::allocate( initSize));
                end_ = begin_ + initSize;
                cur_ = begin_;
            }
        }

        template<typename Allocator>
        CustomMemoryBuffer<Allocator>::CustomMemoryBuffer( const CustomMemoryBuffer &pb )
        {
            if( pb.begin_ != nullptr )
            {
                begin_ = reinterpret_cast<ElementType *>(Allocator::allocate( pb.allocSize()));
                end_ = begin_ + pb.allocSize();
                Allocator::memcpy( begin_, pb.begin_, pb.allocSize());
                cur_ = begin_ + pb.size();
            }
        }

        template<typename Allocator>
        CustomMemoryBuffer<Allocator> &CustomMemoryBuffer<Allocator>::operator=( const CustomMemoryBuffer<Allocator> &pb )
        {
            if( this == &pb )
                return *this;

            CustomMemoryBuffer<Allocator> copy( pb );
            std::swap( cur_, copy.cur_ );
            std::swap( begin_, copy.begin_ );
            std::swap( end_, copy.end_ );
            return *this;
        }

        template<typename Allocator>
        CustomMemoryBuffer<Allocator>::~CustomMemoryBuffer()
        {
            if( begin_ != nullptr )
                Allocator::deallocate( begin_ );
        }

        template<typename Allocator>
        void CustomMemoryBuffer<Allocator>::resize( std::size_t newSize )
        {
            if( newSize > allocSize())
            {
                auto offset = cur_ - begin_;

                ElementType *newBegin;

                newBegin = reinterpret_cast<ElementType *>(Allocator::allocate( newSize ));

                // memcpy: If either dest or src is an invalid or null pointer, the behavior is undefined, even if count is zero.
                if(begin_) {
                    Allocator::memcpy( newBegin, begin_, size_t(end_ - begin_) );
                }

                std::swap( begin_, newBegin );
                if( newBegin != nullptr )
                    Allocator::deallocate( newBegin );

                end_ = begin_ + newSize;
                cur_ = begin_ + offset;
            }

        }

        template<typename Allocator>
        typename CustomMemoryBuffer<Allocator>::ElementType *CustomMemoryBuffer<Allocator>::advance( std::size_t bytes )
        {
            resize( size() + bytes );
            auto result = cur_;
            cur_ += bytes;
            WALBERLA_ASSERT_LESS_EQUAL( cur_, end_ )
            return result;
        }

        template<typename Allocator>
        typename CustomMemoryBuffer<Allocator>::ElementType *CustomMemoryBuffer<Allocator>::advanceNoResize( std::size_t bytes )
        {
            auto newSize = size() + bytes;
            if( newSize <= allocSize())
                return advance( bytes );
            else
                return nullptr;
        }

    }

}

#endif //TURBINECORE_CUSTOMMEMORYBUFFER_H
