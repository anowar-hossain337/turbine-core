
#pragma once

#ifndef TURBINECORE_VECTOR2_H
#define TURBINECORE_VECTOR2_H

#include <cmath>
#include <cassert>
#include <iostream>

#include <core/math/Vector2.h>

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/ProjectDefines.h"

namespace turbine_core {

    namespace math {

        template< typename T >
        class Vector2 {

        public:

            HOST_DEVICE_PREFIX          Vector2() {}
            HOST_DEVICE_PREFIX explicit Vector2( const T & init ) { v_[0] = v_[1] = init; }
            HOST_DEVICE_PREFIX          Vector2( const T & x, const T & y ) { v_[0] = x; v_[1] = y; }
            HOST_DEVICE_PREFIX          Vector2( const Vector2 & v ) { v_[0] = v[0]; v_[1] = v[1]; }
            HOST_DEVICE_PREFIX explicit Vector2( const T * v ) { v_[0] = v[0]; v_[1] = v[1]; }

            HOST_DEVICE_PREFIX explicit Vector2( const walberla::Vector2<T> & vector )
            : Vector2(vector.data())
            {}

            HOST_DEVICE_PREFIX Vector2& operator=(const Vector2 & v) {
                if(this != &v) {v_[0] = v.v_[0]; v_[1] = v.v_[1];}
                return *this;
            }

            HOST_DEVICE_PREFIX bool operator==(const Vector2 & v) const {
                if( fabs(v_[0]-v.v_[0])<1e-8 && fabs(v_[1]-v.v_[1])<1e-8 )
                    return true;
                return false;
            }
            HOST_DEVICE_PREFIX bool operator!=(const Vector2 & v) const { return !(this->operator==(v)); }

            HOST_DEVICE_PREFIX       T & operator[](uint_t idx) { assert(idx < 2); return v_[idx]; }
            HOST_DEVICE_PREFIX const T & operator[](uint_t idx) const  { assert(idx < 2); return v_[idx]; }

            HOST_DEVICE_PREFIX inline Vector2  operator-() const { return Vector2(-v_[0], -v_[1]); }

            HOST_DEVICE_PREFIX inline Vector2& operator+=( const Vector2& rhs ) {
                v_[0] += rhs.v_[0];
                v_[1] += rhs.v_[1];
                return *this;
            }

            HOST_DEVICE_PREFIX inline Vector2& operator-=( const Vector2& rhs ) {
                v_[0] -= rhs.v_[0];
                v_[1] -= rhs.v_[1];
                return *this;
            }

            HOST_DEVICE_PREFIX inline Vector2& operator*=( T rhs ) {
                v_[0] *= rhs; v_[1] *= rhs; return *this;
            }

            HOST_DEVICE_PREFIX inline Vector2& operator/=( T rhs ) {
                v_[0] /= rhs; v_[1] /= rhs; return *this;
            }

            HOST_DEVICE_PREFIX inline Vector2  operator+ ( const Vector2& rhs ) const {
                return Vector2( v_[0] + rhs.v_[0], v_[1] + rhs.v_[1] );
            }

            HOST_DEVICE_PREFIX inline Vector2  operator- ( const Vector2& rhs ) const {
                return Vector2( v_[0] - rhs.v_[0], v_[1] - rhs.v_[1] );
            }

            HOST_DEVICE_PREFIX inline Vector2  operator* ( T rhs ) const {
                return Vector2( v_[0] * rhs, v_[1] * rhs );
            }

            HOST_DEVICE_PREFIX inline T operator* ( const Vector2& rhs ) const {
                return v_[0] * rhs.v_[0] + v_[1] * rhs.v_[1];
            }

            HOST_DEVICE_PREFIX inline Vector2  operator/ ( T rhs ) const {
                return Vector2( v_[0] / rhs, v_[1] / rhs );
            }

            HOST_DEVICE_PREFIX inline T min() const {
                return min(v_[0], v_[1]);
            }

            HOST_DEVICE_PREFIX inline T max() const {
                return max(v_[0], v_[1]);
            }

            HOST_DEVICE_PREFIX inline void set( const T x, const T y ) {
                v_[0] = x;
                v_[1] = y;
            }

            HOST_DEVICE_PREFIX T length() const {
                return sqrt( v_[0] * v_[0] + v_[1] * v_[1] );
            }
            HOST_DEVICE_PREFIX T sqrLength() const  {
                return v_[0] * v_[0] + v_[1] * v_[1];
            }

            HOST_DEVICE_PREFIX void reset() {
                v_[0] = v_[1] = T(0);
            }

            HOST_DEVICE_PREFIX T* data() { return v_; }
            HOST_DEVICE_PREFIX T const * data() const {return v_;}

        private:

            T v_[2]{T(), T()};

        };

        template< typename T >
        HOST_DEVICE_PREFIX inline Vector2<T> normalize( const Vector2<T> & v ) {
            return v / v.length();
        }

        template< typename Type >
        HOST_DEVICE_PREFIX inline bool isnan( const Vector2<Type>& v ) {
            return (v[0] != v[0]) && (v[1] != v[1]);
        }

        template< typename Type >
        HOST_DEVICE_PREFIX inline bool isinf( const Vector2<Type>& v ) {
            return isinf(v[0]) && isinf(v[1]);
        }

        template< typename Type >
        HOST_DEVICE_PREFIX inline bool finite( const Vector2<Type>& v ) {
            return isfinite(v[0]) && isfinite(v[1]);
        }

        template< typename Type >
        HOST_DEVICE_PREFIX inline Vector2<Type> abs( const Vector2<Type>& v ) {
            return Vector2<Type>( abs(v[0]), abs(v[1]) );
        }

        template< typename Type >
        HOST_DEVICE_PREFIX inline Vector2<Type> fabs( const Vector2<Type>& v ) {
            return Vector2<Type>( fabs(v[0]), fabs(v[1]) );
        }

        template< typename Type >
        std::ostream& operator<<( std::ostream& os, const Vector2<Type>& v )
        {
            return os << "<" << v[0] << "," << v[1]  << ">";
        }

        //TODO assertion for valid size
        template< typename Buffer_T, typename Type_T >
        HOST_DEVICE_PREFIX Buffer_T * operator<<( Buffer_T * buffer, const Vector2<Type_T>& v ) {
            constexpr auto size = sizeof(Vector2<Type_T>);
            memcpy(buffer, v.data(), size);
            return buffer + size;
        }

        //TODO assertion for valid size
        template< typename Buffer_T, typename Type_T >
        HOST_DEVICE_PREFIX Buffer_T * operator>>( Buffer_T * buffer, Vector2<Type_T>& v ) {
            constexpr auto size = sizeof(Vector2<Type_T>);
            memcpy(v.data(), buffer, size);
            return buffer + size;
        }

    } // namespace math

    using math::Vector2;

}

#endif //TURBINECORE_VECTOR3_H
