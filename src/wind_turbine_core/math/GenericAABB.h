
#pragma once

#ifndef TURBINECORE_GENERICAABB_H
#define TURBINECORE_GENERICAABB_H

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/math/Overloads.h"

#include "wind_turbine_core/math/Vector3.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include <core/math/AABB.h>

namespace turbine_core {

    namespace math {

        template< typename Type_T >
        class GenericAABB {

        public:

            typedef Type_T              value_type;
            typedef Vector3< Type_T >   vector_type;

            HOST_DEVICE_PREFIX GenericAABB() {}

            HOST_DEVICE_PREFIX GenericAABB( const GenericAABB & other )
                    : minCorner_(other.minCorner_), maxCorner_(other.maxCorner_)
            {}

            HOST_DEVICE_PREFIX GenericAABB( const walberla::AABB & other )
                    : minCorner_(other.minCorner().data()), maxCorner_(other.maxCorner().data())
            {}

            HOST_DEVICE_PREFIX GenericAABB( const vector_type & corner0, const vector_type & corner1 )
            {
                for( uint_t i = 0; i < 3; ++i )
                {
                    if( corner0[i] < corner1[i] )
                    {
                        minCorner_[i] = corner0[i];
                        maxCorner_[i] = corner1[i];
                    }
                    else
                    {
                        minCorner_[i] = corner1[i];
                        maxCorner_[i] = corner0[i];
                    }
                }
            }

            HOST_DEVICE_PREFIX GenericAABB( const value_type x0, const value_type y0, const value_type z0,
                                            const value_type x1, const value_type y1, const value_type z1 ) {
                if( x0 < x1 ) {
                    minCorner_[0] = x0;
                    maxCorner_[0] = x1;
                } else {
                    minCorner_[0] = x1;
                    maxCorner_[0] = x0;
                }

                if( y0 < y1 ) {
                    minCorner_[1] = y0;
                    maxCorner_[1] = y1;
                } else {
                    minCorner_[1] = y1;
                    maxCorner_[1] = y0;
                }

                if( z0 < z1 ) {
                    minCorner_[2] = z0;
                    maxCorner_[2] = z1;
                } else {
                    minCorner_[2] = z1;
                    maxCorner_[2] = z0;
                }
            }

            HOST_DEVICE_PREFIX GenericAABB( const GenericAABB & lhs, const GenericAABB & rhs )
                    : minCorner_( math::max( lhs.minCorner_[0], rhs.minCorner_[0] ),
                                  math::max( lhs.minCorner_[1], rhs.minCorner_[1] ),
                                  math::max( lhs.minCorner_[2], rhs.minCorner_[2] ) ),

                      maxCorner_( math::max( minCorner_[0], math::min( lhs.maxCorner_[0], rhs.maxCorner_[0] ) ),
                                  math::max( minCorner_[1], math::min( lhs.maxCorner_[1], rhs.maxCorner_[1] ) ),
                                  math::max( minCorner_[2], math::min( lhs.maxCorner_[2], rhs.maxCorner_[2] ) ) )
            {}

            template< typename InputIterator >
            HOST_DEVICE_PREFIX GenericAABB( InputIterator first, InputIterator last ) {
                if( first != last ) {
                    merge( first, last );
                } else {
                    init();
                }
            }

            template< typename AABB_T >
            HOST_DEVICE_PREFIX explicit GenericAABB( const AABB_T & aabb )
            : minCorner_(aabb.minCorner().data()), maxCorner_(aabb.maxCorner().data())
            {}

            HOST_DEVICE_PREFIX static GenericAABB createFromMinMaxCorner( const vector_type & theMinCorner, const vector_type & theMaxCorner ) {
                return GenericAABB(theMinCorner, theMaxCorner);
            }

            HOST_DEVICE_PREFIX static GenericAABB createFromMinMaxCorner( const value_type minX, const value_type minY, const value_type minZ,
                                                                          const value_type maxX, const value_type maxY, const value_type maxZ ) {
                return GenericAABB( minX, minY, minZ,
                                    maxX, maxY, maxZ );
            }

            HOST_DEVICE_PREFIX GenericAABB & operator=( const GenericAABB & other ) {
                minCorner_ = other.minCorner_;
                maxCorner_ = other.maxCorner_;
                return *this;
            }

            HOST_DEVICE_PREFIX const vector_type & minCorner() const {
                return minCorner_;
            }

            HOST_DEVICE_PREFIX const vector_type & maxCorner() const {
                return maxCorner_;
            }

            HOST_DEVICE_PREFIX const vector_type & min() const { return minCorner(); }
            HOST_DEVICE_PREFIX const vector_type & max() const { return maxCorner(); }

            HOST_DEVICE_PREFIX value_type xMin() const {
                return minCorner_[0];
            }

            HOST_DEVICE_PREFIX value_type yMin() const {
                return minCorner_[1];
            }

            HOST_DEVICE_PREFIX value_type zMin() const {
                return minCorner_[2];
            }

            HOST_DEVICE_PREFIX value_type min( const uint_t index ) const {
                assert( index < uint_t(3) );
                return minCorner_[index];
            }

            HOST_DEVICE_PREFIX value_type xMax() const {
                return maxCorner_[0];
            }

            HOST_DEVICE_PREFIX value_type yMax() const {
                return maxCorner_[1];
            }

            HOST_DEVICE_PREFIX value_type zMax() const {
                return maxCorner_[2];
            }

            HOST_DEVICE_PREFIX value_type max( const uint_t index ) const {
                assert( index < uint_t(3) );
                return maxCorner_[index];
            }

            HOST_DEVICE_PREFIX bool empty() const {
                return minCorner_ == maxCorner_;
            }

            HOST_DEVICE_PREFIX vector_type sizes() const {
                return maxCorner_ - minCorner_;
            }

            HOST_DEVICE_PREFIX value_type size( const uint_t i ) const {
                assert( i < uint_t(3) );
                return maxCorner_[i] - minCorner_[i];
            }

            HOST_DEVICE_PREFIX value_type xSize() const {
                return size(0);
            }

            HOST_DEVICE_PREFIX value_type ySize() const {
                return size(1);
            }

            HOST_DEVICE_PREFIX value_type zSize() const {
                return size(2);
            }

            HOST_DEVICE_PREFIX value_type volume() const {
                return xSize() * ySize() * zSize();
            }

            HOST_DEVICE_PREFIX vector_type center() const {
                return (minCorner_ + maxCorner_) * value_type(0.5);
            }

            HOST_DEVICE_PREFIX bool contains( const value_type x, const value_type y, const value_type z ) const {
                return x >= minCorner_[0] && x < maxCorner_[0] &&
                       y >= minCorner_[1] && y < maxCorner_[1] &&
                       z >= minCorner_[2] && z < maxCorner_[2];
            }

            HOST_DEVICE_PREFIX bool contains( const vector_type & point ) const {
                return contains(point[0], point[1], point[2]);
            }

            template<typename Other>
            HOST_DEVICE_PREFIX bool contains( const Vector3<Other> & point ) const {
                return contains((Type_T)point[0], (Type_T)point[1], (Type_T)point[2]);
            }

            HOST_DEVICE_PREFIX bool contains( const value_type x, const value_type y, const value_type z, const value_type dx ) const
            {
                return x >= ( minCorner_[0] - dx ) && x < ( maxCorner_[0] + dx ) &&
                       y >= ( minCorner_[1] - dx ) && y < ( maxCorner_[1] + dx ) &&
                       z >= ( minCorner_[2] - dx ) && z < ( maxCorner_[2] + dx );
            }

            HOST_DEVICE_PREFIX bool contains( const vector_type & point, const value_type dx ) {
                return contains( point[0], point[1], point[2], dx );
            }

            template<typename Other>
            HOST_DEVICE_PREFIX bool contains( const Vector3<Other> & point, const value_type dx ) {
                return contains( (Type_T)point[0], (Type_T)point[1], (Type_T)point[2], dx );
            }

            HOST_DEVICE_PREFIX bool contains( const GenericAABB & other ) const
            {
                return other.minCorner_[0] >= minCorner_[0] && other.maxCorner_[0] <= maxCorner_[0] &&
                       other.minCorner_[1] >= minCorner_[1] && other.maxCorner_[1] <= maxCorner_[1] &&
                       other.minCorner_[2] >= minCorner_[2] && other.maxCorner_[2] <= maxCorner_[2];
            }

            HOST_DEVICE_PREFIX bool contains( const walberla::AABB & other ) const
            {
                const auto minOther = other.minCorner();
                const auto maxOther = other.maxCorner();

                return minOther[0] >= minCorner_[0] && maxOther[0] <= maxCorner_[0] &&
                       minOther[1] >= minCorner_[1] && maxOther[1] <= maxCorner_[1] &&
                       minOther[2] >= minCorner_[2] && maxOther[2] <= maxCorner_[2];
            }

            HOST_DEVICE_PREFIX bool contains( const GenericAABB & other, const value_type dx ) const
            {
                return other.minCorner_[0] >= ( minCorner_[0] - dx ) && other.maxCorner_[0] <= ( maxCorner_[0] + dx ) &&
                       other.minCorner_[1] >= ( minCorner_[1] - dx ) && other.maxCorner_[1] <= ( maxCorner_[1] + dx ) &&
                       other.minCorner_[2] >= ( minCorner_[2] - dx ) && other.maxCorner_[2] <= ( maxCorner_[2] + dx );
            }

            HOST_DEVICE_PREFIX bool contains( const walberla::AABB & other, const value_type dx  ) const
            {
                const auto minOther = other.minCorner();
                const auto maxOther = other.maxCorner();

                return minOther[0] >= ( minCorner_[0] - dx ) && maxOther[0] <= ( maxCorner_[0] + dx ) &&
                       minOther[1] >= ( minCorner_[1] - dx ) && maxOther[1] <= ( maxCorner_[1] + dx ) &&
                       minOther[2] >= ( minCorner_[2] - dx ) && maxOther[2] <= ( maxCorner_[2] + dx );
            }

            HOST_DEVICE_PREFIX bool containsClosedInterval( const vector_type & point ) const
            {
                return point[0] >= minCorner_[0] && point[0] <= maxCorner_[0] &&
                       point[1] >= minCorner_[1] && point[1] <= maxCorner_[1] &&
                       point[2] >= minCorner_[2] && point[2] <= maxCorner_[2];
            }

            HOST_DEVICE_PREFIX bool containsClosedInterval( const vector_type & point, const value_type dx ) const
            {
                return point[0] >= ( minCorner_[0] - dx ) && point[0] <= ( maxCorner_[0] + dx ) &&
                       point[1] >= ( minCorner_[1] - dx ) && point[1] <= ( maxCorner_[1] + dx ) &&
                       point[2] >= ( minCorner_[2] - dx ) && point[2] <= ( maxCorner_[2] + dx );
            }

            HOST_DEVICE_PREFIX GenericAABB getExtended( const value_type eps ) const
            {
                vector_type newMinCorner( minCorner_[0] - eps, minCorner_[1] - eps, minCorner_[2] - eps );

                return createFromMinMaxCorner( newMinCorner[0], newMinCorner[1], newMinCorner[2],
                                               math::max( newMinCorner[0], maxCorner_[0] + eps ),
                                               math::max( newMinCorner[1], maxCorner_[1] + eps ),
                                               math::max( newMinCorner[2], maxCorner_[2] + eps ) );
            }

            HOST_DEVICE_PREFIX GenericAABB getExtended( const vector_type & eps ) const
            {
                vector_type newMinCorner( minCorner_ - eps );
                return createFromMinMaxCorner( newMinCorner[0], newMinCorner[1], newMinCorner[2],
                                               math::max( newMinCorner[0], maxCorner_[0] + eps[0] ),
                                               math::max( newMinCorner[1], maxCorner_[1] + eps[1] ),
                                               math::max( newMinCorner[2], maxCorner_[2] + eps[2] ) );
            }

            HOST_DEVICE_PREFIX GenericAABB getTranslated( const vector_type & translation ) const
            {
                return createFromMinMaxCorner( minCorner_ + translation, maxCorner_ + translation );
            }

            HOST_DEVICE_PREFIX GenericAABB getScaled( const value_type factor ) const
            {
                vector_type theCenter = center();

                GenericAABB result = getTranslated( -theCenter );
                result.minCorner_ *= factor;
                result.maxCorner_ *= factor;
                result.translate( theCenter );

                return result;
            }

            HOST_DEVICE_PREFIX GenericAABB getScaled( const vector_type & factors ) const
            {
                vector_type theCenter = center();

                GenericAABB result = getTranslated( -theCenter );

                for( uint_t i = 0; i < 3; ++i ) {
                    result.minCorner_[i] *= factors[i];
                    result.maxCorner_[i] *= factors[i];
                }

                result.translate( theCenter );

                return result;
            }

            HOST_DEVICE_PREFIX GenericAABB getMerged( const vector_type & point ) const
            {
                return createFromMinMaxCorner( math::min( minCorner_[0], point[0] ),
                                               math::min( minCorner_[1], point[1] ),
                                               math::min( minCorner_[2], point[2] ),
                                               math::max( maxCorner_[0], point[0] ),
                                               math::max( maxCorner_[1], point[1] ),
                                               math::max( maxCorner_[2], point[2] ) );
            }

            HOST_DEVICE_PREFIX GenericAABB getMerged( const GenericAABB & other ) const
            {
                return createFromMinMaxCorner( math::min( minCorner_[0], other.minCorner_[0] ),
                                               math::min( minCorner_[1], other.minCorner_[1] ),
                                               math::min( minCorner_[2], other.minCorner_[2] ),
                                               math::max( maxCorner_[0], other.maxCorner_[0] ),
                                               math::max( maxCorner_[1], other.maxCorner_[1] ),
                                               math::max( maxCorner_[2], other.maxCorner_[2] ) );
            }

            HOST_DEVICE_PREFIX GenericAABB getMerged( const walberla::AABB & other ) const
            {
                const auto minOther = other.minCorner();
                const auto maxOther = other.maxCorner();

                return createFromMinMaxCorner( math::min( minCorner_[0], minOther[0] ),
                                               math::min( minCorner_[1], minOther[1] ),
                                               math::min( minCorner_[2], minOther[2] ),
                                               math::max( maxCorner_[0], maxOther[0] ),
                                               math::max( maxCorner_[1], maxOther[1] ),
                                               math::max( maxCorner_[2], maxOther[2] ) );
            }

            HOST_DEVICE_PREFIX bool intersects( const GenericAABB & other ) const
            {
                return other.maxCorner_[0] > minCorner_[0] && other.minCorner_[0] < maxCorner_[0] &&
                       other.maxCorner_[1] > minCorner_[1] && other.minCorner_[1] < maxCorner_[1] &&
                       other.maxCorner_[2] > minCorner_[2] && other.minCorner_[2] < maxCorner_[2] &&
                       !empty() && !other.empty();
            }

            HOST_DEVICE_PREFIX bool intersects( const walberla::AABB & other ) const
            {
                const auto minOther = other.minCorner();
                const auto maxOther = other.maxCorner();

                return maxOther[0] > minCorner_[0] && minOther[0] < maxCorner_[0] &&
                       maxOther[1] > minCorner_[1] && minOther[1] < maxCorner_[1] &&
                       maxOther[2] > minCorner_[2] && minOther[2] < maxCorner_[2] &&
                       !empty() && !other.empty();
            }

            HOST_DEVICE_PREFIX bool intersects( const GenericAABB & other, const value_type dx ) const
            {
                assert( minCorner_[0] - dx < maxCorner_[0] + dx );
                assert( minCorner_[1] - dx < maxCorner_[1] + dx );
                assert( minCorner_[2] - dx < maxCorner_[2] + dx );

                return other.maxCorner_[0] > ( minCorner_[0] - dx ) && other.minCorner_[0] < ( maxCorner_[0] + dx ) &&
                       other.maxCorner_[1] > ( minCorner_[1] - dx ) && other.minCorner_[1] < ( maxCorner_[1] + dx ) &&
                       other.maxCorner_[2] > ( minCorner_[2] - dx ) && other.minCorner_[2] < ( maxCorner_[2] + dx ) &&
                       !other.empty() &&
                       minCorner_[0] - dx < maxCorner_[0] + dx &&
                       minCorner_[1] - dx < maxCorner_[1] + dx &&
                       minCorner_[2] - dx < maxCorner_[2] + dx;
            }

            HOST_DEVICE_PREFIX bool intersects( const walberla::AABB & other, const value_type dx ) const
            {
                const auto minOther = other.minCorner();
                const auto maxOther = other.maxCorner();

                assert( minCorner_[0] - dx < maxCorner_[0] + dx );
                assert( minCorner_[1] - dx < maxCorner_[1] + dx );
                assert( minCorner_[2] - dx < maxCorner_[2] + dx );

                return maxOther[0] > ( minCorner_[0] - dx ) && minOther[0] < ( maxCorner_[0] + dx ) &&
                       maxOther[1] > ( minCorner_[1] - dx ) && minOther[1] < ( maxCorner_[1] + dx ) &&
                       maxOther[2] > ( minCorner_[2] - dx ) && minOther[2] < ( maxCorner_[2] + dx ) &&
                       !other.empty() &&
                       minCorner_[0] - dx < maxCorner_[0] + dx &&
                       minCorner_[1] - dx < maxCorner_[1] + dx &&
                       minCorner_[2] - dx < maxCorner_[2] + dx;
            }

            HOST_DEVICE_PREFIX bool intersectsClosedInterval( const GenericAABB & other ) const
            {
                return other.maxCorner_[0] >= minCorner_[0] && other.minCorner_[0] <= maxCorner_[0] &&
                       other.maxCorner_[1] >= minCorner_[1] && other.minCorner_[1] <= maxCorner_[1] &&
                       other.maxCorner_[2] >= minCorner_[2] && other.minCorner_[2] <= maxCorner_[2];
            }

            HOST_DEVICE_PREFIX bool intersectsClosedInterval( const walberla::AABB & other ) const
            {
                const auto minOther = other.minCorner();
                const auto maxOther = other.maxCorner();

                return maxOther[0] >= minCorner_[0] && minOther[0] <= maxCorner_[0] &&
                       maxOther[1] >= minCorner_[1] && minOther[1] <= maxCorner_[1] &&
                       maxOther[2] >= minCorner_[2] && minOther[2] <= maxCorner_[2];
            }

            HOST_DEVICE_PREFIX bool intersectsClosedInterval( const GenericAABB & other, const value_type dx ) const
            {
                return other.maxCorner_[0] >= ( minCorner_[0] - dx ) && other.minCorner_[0] <= ( maxCorner_[0] + dx ) &&
                       other.maxCorner_[1] >= ( minCorner_[1] - dx ) && other.minCorner_[1] <= ( maxCorner_[1] + dx ) &&
                       other.maxCorner_[2] >= ( minCorner_[2] - dx ) && other.minCorner_[2] <= ( maxCorner_[2] + dx );
            }

            HOST_DEVICE_PREFIX bool intersectsClosedInterval( const walberla::AABB & other, const value_type dx ) const
            {
                const auto minOther = other.minCorner();
                const auto maxOther = other.maxCorner();

                return maxOther[0] >= ( minCorner_[0] - dx ) && minOther[0] <= ( maxCorner_[0] + dx ) &&
                       maxOther[1] >= ( minCorner_[1] - dx ) && minOther[1] <= ( maxCorner_[1] + dx ) &&
                       maxOther[2] >= ( minCorner_[2] - dx ) && minOther[2] <= ( maxCorner_[2] + dx );
            }

            value_type intersectionVolume( const GenericAABB & other ) const
            {
                return getIntersection( other ).volume();
            }

            value_type intersectionVolume( const walberla::AABB & other ) const
            {
                return getIntersection( other ).volume();
            }

            HOST_DEVICE_PREFIX GenericAABB getIntersection( const GenericAABB & other ) const
            {
                return GenericAABB( *this, other );
            }

            HOST_DEVICE_PREFIX GenericAABB getIntersection( const walberla::AABB & other ) const
            {
                return GenericAABB( *this, other );
            }

            HOST_DEVICE_PREFIX bool isIdentical( const GenericAABB & other ) const
            {
                return (minCorner_ == other.minCorner_) && (maxCorner_ == other.maxCorner_);
            }

            HOST_DEVICE_PREFIX bool isIdentical( const walberla::AABB & other ) const
            {
                return (minCorner_ == other.minCorner()) && (maxCorner_ == other.maxCorner());
            }

            HOST_DEVICE_PREFIX bool isEqual( const GenericAABB & other ) const {
                return (minCorner_ == other.minCorner_) && (maxCorner_ == other.maxCorner_);
            }

            HOST_DEVICE_PREFIX bool isEqual( const walberla::AABB & other ) const {
                return (minCorner_ == other.minCorner()) && (maxCorner_ == other.maxCorner());
            }

            HOST_DEVICE_PREFIX value_type sqDistance( const vector_type & point ) const
            {
                const vector_type d = point - minCorner_;
                const vector_type theSizes = sizes();
                value_type sqDist(0);

                for( uint_t i = 0; i < 3; i++ ) {
                    if( d[i] < 0 ) {
                        sqDist += d[i] * d[i];
                    } else if( d[i] > theSizes[i] ) {
                        const value_type axisDist = d[i] - theSizes[i];
                        sqDist += axisDist * axisDist;
                    }
                }

                assert( sqDist >= value_type( 0 ) );

                return sqDist;
            }

            HOST_DEVICE_PREFIX value_type sqSignedDistance( const vector_type & point ) const
            {
                const vector_type d        = point - minCorner_;
                const vector_type theSizes = sizes();

                bool inside = d[0] >= value_type( 0 ) && d[0] < theSizes[0] &&
                              d[1] >= value_type( 0 ) && d[1] < theSizes[1] &&
                              d[2] >= value_type( 0 ) && d[2] < theSizes[2];

                if( !inside )
                    return sqDistance( point );

                value_type sqAxisDist[3];

                for( uint_t i = 0; i < 3; ++i ) {
                    assert( d[i] >= 0 );
                    assert( d[i] < theSizes[i] );

                    if( d[i] < theSizes[i] * value_type(0.5) ) {
                        sqAxisDist[i] = d[i] * d[i];
                    } else {
                        const value_type axisDist = d[i] - theSizes[i];
                        sqAxisDist[i] = axisDist * axisDist;
                    }
                }

                return -math::min( sqAxisDist[0], math::min( sqAxisDist[1], sqAxisDist[2] ) );
            }

            HOST_DEVICE_PREFIX value_type sqMaxDistance( const vector_type & point ) const
            {
                const vector_type d0 = point - minCorner_;
                const vector_type d1 = point - maxCorner_;
                value_type sqDist(0);

                for( uint_t i = 0; i < 3; i++ )
                {
                    if( fabs( d0[i] ) > fabs( d1[i] ) )
                        sqDist += d0[i] * d0[i];
                    else
                        sqDist += d1[i] * d1[i];
                }

                WALBERLA_ASSERT_GREATER_EQUAL( sqDist, value_type( 0 ) );

                return sqDist;
            }

            HOST_DEVICE_PREFIX value_type distance( const vector_type & point ) const
            {
                return sqrt( sqDistance( point ) );
            }

            HOST_DEVICE_PREFIX value_type signedDistance( const vector_type & point ) const
            {
                const vector_type d        = point - minCorner_;
                const vector_type theSizes = sizes();

                bool inside = d[0] >= value_type( 0 ) && d[0] < theSizes[0] &&
                              d[1] >= value_type( 0 ) && d[1] < theSizes[1] &&
                              d[2] >= value_type( 0 ) && d[2] < theSizes[2];

                if( !inside )
                    return distance( point );

                value_type axisDist[3];

                for( uint_t i = 0; i < 3; ++i )
                {
                    WALBERLA_ASSERT_GREATER_EQUAL( d[i], 0 );
                    WALBERLA_ASSERT_LESS( d[i], theSizes[i] );

                    if( d[i] < theSizes[i] * value_type( 0.5 ) )
                    {
                        axisDist[i] = d[i];
                    }
                    else
                    {
                        axisDist[i] = theSizes[i] - d[i];
                    }
                }

                return -math::min( axisDist[0], math::min( axisDist[1], axisDist[2] ) );
            }

            HOST_DEVICE_PREFIX value_type maxDistance( const vector_type & point ) const
            {
                return sqrt( sqMaxDistance( point ) );
            }

            HOST_DEVICE_PREFIX value_type sqDistance( const GenericAABB & other ) const
            {
                value_type theSqDistance = value_type(0);

                if( other.maxCorner_[0] < minCorner_[0] ) {
                    theSqDistance += sq( minCorner_[0] - other.maxCorner_[0] );
                } else if( other.minCorner_[0] > maxCorner_[0] ) {
                    theSqDistance += sq( other.minCorner_[0] - maxCorner_[0] );
                }

                if( other.maxCorner_[1] < minCorner_[1] ) {
                    theSqDistance += sq( minCorner_[1] - other.maxCorner_[1] );
                } else if( other.minCorner_[1] > maxCorner_[1] ) {
                    theSqDistance += sq( other.minCorner_[1] - maxCorner_[1] );
                }

                if( other.maxCorner_[2] < minCorner_[2] ) {
                    theSqDistance += sq( minCorner_[2] - other.maxCorner_[2] );
                } else if( other.minCorner_[2] > maxCorner_[2] ) {
                    theSqDistance += sq( other.minCorner_[2] - maxCorner_[2] );
                }

                assert( theSqDistance >= value_type(0) );
                assert( !intersects( other ) || ( theSqDistance == value_type(0) ) ); // intersect => distance == 0

                return theSqDistance;
            }

            HOST_DEVICE_PREFIX value_type sqMaxDistance( const GenericAABB & other ) const
            {
                value_type theSqMaxDistance = value_type(0);

                theSqMaxDistance += sq( max( maxCorner_[0] - other.minCorner_[0], other.maxCorner_[0] - minCorner_[0] ) );
                theSqMaxDistance += sq( max( maxCorner_[1] - other.minCorner_[1], other.maxCorner_[1] - minCorner_[1] ) );
                theSqMaxDistance += sq( max( maxCorner_[2] - other.minCorner_[2], other.maxCorner_[2] - minCorner_[2] ) );

                assert( theSqMaxDistance >= value_type(0) );

                return theSqMaxDistance;
            }

            HOST_DEVICE_PREFIX void init()
            {
                minCorner_.set( value_type( 0 ), value_type( 0 ), value_type( 0 ) );
                maxCorner_.set( value_type( 0 ), value_type( 0 ), value_type( 0 ) );

                assert( empty() );
            }

            HOST_DEVICE_PREFIX void init( const vector_type & corner0, const vector_type & corner1 )
            {
                for( uint_t i = 0; i < 3; ++i )
                {
                    if( corner0[i] < corner1[i] )
                    {
                        minCorner_[i] = corner0[i];
                        maxCorner_[i] = corner1[i];
                    }
                    else
                    {
                        minCorner_[i] = corner1[i];
                        maxCorner_[i] = corner0[i];
                    }
                }

            }

            HOST_DEVICE_PREFIX void init( const value_type x0, const value_type y0, const value_type z0,
                                          const value_type x1, const value_type y1, const value_type z1 )
            {
                if( x0 < x1 ) {
                    minCorner_[0] = x0;
                    maxCorner_[0] = x1;
                } else {
                    minCorner_[0] = x1;
                    maxCorner_[0] = x0;
                }

                if( y0 < y1 ) {
                    minCorner_[1] = y0;
                    maxCorner_[1] = y1;
                } else {
                    minCorner_[1] = y1;
                    maxCorner_[1] = y0;
                }

                if( z0 < z1 ) {
                    minCorner_[2] = z0;
                    maxCorner_[2] = z1;
                } else {
                    minCorner_[2] = z1;
                    maxCorner_[2] = z0;
                }

            }

            HOST_DEVICE_PREFIX void initMinMaxCorner( const value_type minX, const value_type minY, const value_type minZ,
                                                      const value_type maxX, const value_type maxY, const value_type maxZ )
            {
                minCorner_[0] = minX;
                minCorner_[1] = minY;
                minCorner_[2] = minZ;

                maxCorner_[0] = maxX;
                maxCorner_[1] = maxY;
                maxCorner_[2] = maxZ;
            }

            HOST_DEVICE_PREFIX void initMinMaxCorner( const vector_type & theMinCorner, const vector_type & theMaxCorner )
            {
                minCorner_ = theMinCorner;
                maxCorner_ = theMaxCorner;
            }

            HOST_DEVICE_PREFIX void setAxisBounds( const uint_t index, const value_type value1, const value_type value2 )
            {
                assert( index < 3 );

                if ( value1 < value2 ) {
                    minCorner_[index] = value1;
                    maxCorner_[index] = value2;
                } else {
                    minCorner_[index] = value2;
                    maxCorner_[index] = value1;
                }
            }

            HOST_DEVICE_PREFIX void extend( const value_type eps )
            {
                minCorner_[0] -= eps;
                minCorner_[1] -= eps;
                minCorner_[2] -= eps;

                maxCorner_[0] += eps;
                maxCorner_[1] += eps;
                maxCorner_[2] += eps;

                for( uint_t i = 0; i < 3; ++i ) {
                    if (minCorner_[i] > maxCorner_[i])
                        maxCorner_[i] = minCorner_[i];
                }
            }

            HOST_DEVICE_PREFIX void extend( const vector_type & eps )
            {
                minCorner_ -= eps;
                maxCorner_ += eps;

                for( uint_t i = 0; i < 3; ++i ) {
                    if(minCorner_[i] > maxCorner_[i])
                        maxCorner_[i] = minCorner_[i];
                }
            }

            HOST_DEVICE_PREFIX void setCenter( const vector_type & center )
            {
                translate(center - this->center());
            }

            HOST_DEVICE_PREFIX void translate( const vector_type & translation )
            {
                minCorner_ += translation;
                maxCorner_ += translation;
            }

            HOST_DEVICE_PREFIX void scale( const value_type factor )
            {
                assert( factor > value_type( 0 ) );

                const vector_type theCenter = center();

                translate( -theCenter );

                minCorner_ *= factor;
                maxCorner_ *= factor;

                translate( theCenter );
            }

            HOST_DEVICE_PREFIX void scale( const vector_type & factors )
            {
                assert( factors[0] > value_type( 0 ) );
                assert( factors[1] > value_type( 0 ) );
                assert( factors[2] > value_type( 0 ) );

                const vector_type theCenter = center();

                translate( -theCenter );

                for( uint_t i = 0; i < 3; ++i ) {
                    minCorner_[i] *= factors[i];
                    maxCorner_[i] *= factors[i];
                }

                translate( theCenter );
            }

            HOST_DEVICE_PREFIX void merge( const vector_type & point ) {
                for( uint_t i = 0; i < 3; ++i ) {
                    if( point[i] < minCorner_[i] )
                        minCorner_[i] = point[i];
                    if( point[i] > maxCorner_[i] )
                        maxCorner_[i] = point[i];
                }
            }

            HOST_DEVICE_PREFIX void merge( const GenericAABB & other ) {
                for( uint_t i = 0; i < 3; ++i ) {
                    if( other.minCorner_[i] < minCorner_[i] )
                        minCorner_[i] = other.minCorner_[i];
                    if( other.maxCorner_[i] > maxCorner_[i] )
                        maxCorner_[i] = other.maxCorner_[i];
                }
            }

            template< typename InputIterator >
            HOST_DEVICE_PREFIX void merge( InputIterator first, InputIterator last )
            {
                while( first != last )
                    merge( *first++ );
            }

            HOST_DEVICE_PREFIX void intersect( const GenericAABB & other )
            {
                minCorner_[0] = math::max( minCorner_[0], other.minCorner_[0] );
                minCorner_[1] = math::max( minCorner_[1], other.minCorner_[1] );
                minCorner_[2] = math::max( minCorner_[2], other.minCorner_[2] );

                maxCorner_[0] = math::max( minCorner_[0], math::min( maxCorner_[0], other.maxCorner_[0] ) );
                maxCorner_[1] = math::max( minCorner_[1], math::min( maxCorner_[1], other.maxCorner_[1] ) );
                maxCorner_[2] = math::max( minCorner_[2], math::min( maxCorner_[2], other.maxCorner_[2] ) );
            }

        private:

            vector_type minCorner_{};
            vector_type maxCorner_{};

        };

        template< typename T, typename U >
        HOST_DEVICE_PREFIX bool operator==( const GenericAABB< T > & lhs, const GenericAABB< U > & rhs )
        {
            return lhs.isEqual( rhs );
        }

        template< typename T, typename U >
        HOST_DEVICE_PREFIX bool operator!=( const GenericAABB< T > & lhs, const GenericAABB< U > & rhs )
        {
            return !lhs.isEqual( rhs );
        }

        template< typename T >
        HOST_PREFIX std::ostream& operator<<( std::ostream& os, const GenericAABB< T > & aabb )
        {
            return os << "[ " << aabb.minCorner() << ", " << aabb.maxCorner() << " ]";
        }

        template< typename T >
        HOST_PREFIX std::istream& operator>>( std::istream& is, GenericAABB< T > & aabb )
        {
            if( !is ) return is;

            char bracket0, bracket1;
            char comma;
            typename GenericAABB< T >::vector_type corner0, corner1;

            const std::istream::pos_type pos( is.tellg() );
            const std::istream::fmtflags oldFlags( is.flags() );

            // Setting the 'skip whitespaces' flag
            is >> std::skipws;

            // Extracting the vector
            if( !( is >> bracket0 >> corner0 >> comma >> corner1 >> bracket1 ) ||
                comma != ',' || bracket0 != '[' || bracket1 != ']' )
            {
                is.clear();
                is.seekg( pos );
                is.setstate( std::istream::failbit );
                is.flags( oldFlags );
                return is;
            }

            // Transferring the input to the aabb values
            aabb.init( corner0, corner1 );

            // Resetting the flags
            is.flags( oldFlags );

            return is;
        }

    }

}

#endif //TURBINECORE_GENERICAABB_H
