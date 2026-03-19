
#pragma once

#ifndef TURBINECORE_VECTOR3_H
#define TURBINECORE_VECTOR3_H

#include <cmath>
#include <cassert>
#include <iostream>

#include <core/DataTypes.h>
#include <core/math/Vector3.h>

#include "wind_turbine_core/ProjectDefines.h"

namespace turbine_core {

    namespace math {

        using walberla::real_t;
        using walberla::uint_t;

        template< typename T >
        class Vector3 {

        public:

            HOST_DEVICE_PREFIX constexpr          Vector3() {}
            HOST_DEVICE_PREFIX constexpr explicit Vector3( const T & init ) { v_[0] = v_[1] = v_[2] = init; }
            HOST_DEVICE_PREFIX constexpr          Vector3( const T & x, const T & y, const T & z ) { v_[0] = x; v_[1] = y; v_[2] = z; }
            HOST_DEVICE_PREFIX constexpr          Vector3( const Vector3 & v ) { v_[0] = v[0]; v_[1] = v[1]; v_[2] = v[2]; }
            HOST_DEVICE_PREFIX constexpr explicit Vector3( const T * v ) { v_[0] = v[0]; v_[1] = v[1]; v_[2] = v[2]; }

            HOST_DEVICE_PREFIX Vector3 ( const walberla::Vector3<T> & vector )
            : Vector3(vector.data())
            {}

            HOST_DEVICE_PREFIX Vector3& operator=(const Vector3 & v) {
                if(this != &v) {v_[0] = v.v_[0]; v_[1] = v.v_[1]; v_[2] = v.v_[2];}
                return *this;
            }

            HOST_DEVICE_PREFIX Vector3& operator=(const walberla::Vector3<T> & v) {
                {v_[0] = v[0]; v_[1] = v[1]; v_[2] = v[2];}
                return *this;
            }

            HOST_DEVICE_PREFIX bool operator==(const Vector3 & v) const {
                if( fabs(v_[0]-v.v_[0])<1e-8 && fabs(v_[1]-v.v_[1])<1e-8 && fabs(v_[2]-v.v_[2])<1e-8 )
                    return true;
                return false;
            }
            HOST_DEVICE_PREFIX bool operator!=(const Vector3 & v) const { return !(this->operator==(v)); }

            HOST_DEVICE_PREFIX       T & operator[](uint_t idx) { assert(idx < 3); return v_[idx]; }
            HOST_DEVICE_PREFIX const T & operator[](uint_t idx) const  { assert(idx < 3); return v_[idx]; }

            HOST_DEVICE_PREFIX inline Vector3  operator-() const { return Vector3(-v_[0], -v_[1], -v_[2]); }
            HOST_DEVICE_PREFIX inline Vector3& operator%=( const Vector3& rhs ) {       //cross product
                T tmp0 = v_[1] * rhs.v_[2] - v_[2] * rhs.v_[1];
                T tmp1 = v_[2] * rhs.v_[0] - v_[0] * rhs.v_[2];
                v_[2]  = v_[0] * rhs.v_[1] - v_[1] * rhs.v_[0];
                v_[1]  = tmp1;
                v_[0]  = tmp0;
                return *this;
            }

            HOST_DEVICE_PREFIX inline Vector3& operator+=( const Vector3& rhs ) {
                v_[0] += rhs.v_[0];
                v_[1] += rhs.v_[1];
                v_[2] += rhs.v_[2];
                return *this;
            }

            HOST_DEVICE_PREFIX inline Vector3& operator-=( const Vector3& rhs ) {
                v_[0] -= rhs.v_[0];
                v_[1] -= rhs.v_[1];
                v_[2] -= rhs.v_[2];
                return *this;
            }

            HOST_DEVICE_PREFIX inline Vector3& operator*=( T rhs ) {
                v_[0] *= rhs; v_[1] *= rhs; v_[2] *= rhs; return *this;
            }

            HOST_DEVICE_PREFIX inline Vector3& operator/=( T rhs ) {
                v_[0] /= rhs; v_[1] /= rhs; v_[2] /= rhs; return *this;
            }

            HOST_DEVICE_PREFIX inline Vector3  operator% ( const Vector3& rhs ) const { //cross product
                return Vector3( v_[1] * rhs.v_[2] - v_[2] * rhs.v_[1],
                                v_[2] * rhs.v_[0] - v_[0] * rhs.v_[2],
                                v_[0] * rhs.v_[1] - v_[1] * rhs.v_[0] );
            }

            HOST_DEVICE_PREFIX inline Vector3  operator+ ( const Vector3& rhs ) const {
                return Vector3( v_[0] + rhs.v_[0], v_[1] + rhs.v_[1], v_[2] + rhs.v_[2] );
            }

            HOST_DEVICE_PREFIX inline Vector3  operator- ( const Vector3& rhs ) const {
                return Vector3( v_[0] - rhs.v_[0], v_[1] - rhs.v_[1], v_[2] - rhs.v_[2] );
            }

            HOST_DEVICE_PREFIX inline Vector3  operator* ( T rhs ) const {
                return Vector3( v_[0] * rhs, v_[1] * rhs, v_[2] * rhs );
            }

            HOST_DEVICE_PREFIX inline T operator* ( const Vector3& rhs ) const {
                return v_[0] * rhs.v_[0] + v_[1] * rhs.v_[1] + v_[2] * rhs.v_[2];
            }

            HOST_DEVICE_PREFIX inline Vector3  operator/ ( T rhs ) const {
                return Vector3( v_[0] / rhs, v_[1] / rhs, v_[2] / rhs );
            }

            HOST_DEVICE_PREFIX inline T min() const {
                return min(v_[0], min(v_[1], v_[2]));
            }

            HOST_DEVICE_PREFIX inline T max() const {
                return max(v_[0], max(v_[1], v_[2]));
            }

            HOST_DEVICE_PREFIX inline void set( const T x, const T y, const T z )
            {
                v_[0] = x;
                v_[1] = y;
                v_[2] = z;
            }

            HOST_DEVICE_PREFIX T length() const {
                return sqrt( v_[0] * v_[0] + v_[1] * v_[1] + v_[2] * v_[2] );
            }
            HOST_DEVICE_PREFIX T sqrLength() const  {
                return v_[0] * v_[0] + v_[1] * v_[1] + v_[2] * v_[2];
            }

            HOST_DEVICE_PREFIX void reset() {
                v_[0] = v_[1] = v_[2] = T(0);
            }

            HOST_DEVICE_PREFIX T* data() { return v_; }
            HOST_DEVICE_PREFIX T const * data() const {return v_;}

        private:

            T v_[3]{T(), T(), T()};

        };

        template<typename T>
        HOST_DEVICE_PREFIX inline Vector3<T> operator*(const T t, const Vector3<T> & v) {
            return v * t;
        }

        template< typename T >
        HOST_DEVICE_PREFIX inline Vector3<T> normalize( const Vector3<T> & v ) {
            return v / v.length();
        }

        template< typename Type >
        HOST_DEVICE_PREFIX inline bool isnan( const Vector3<Type>& v ) {
            return (v[0]!=v[0]) && (v[1]!=v[1]) && (v[2]!=v[2]);
        }

        template< typename Type >
        HOST_DEVICE_PREFIX inline bool isinf( const Vector3<Type>& v ) {
            return isinf(v[0]) && isinf(v[1]) && isinf(v[2]);
        }

        template< typename Type >
        HOST_DEVICE_PREFIX inline bool finite( const Vector3<Type>& v ) {
            return isfinite(v[0]) && isfinite(v[1]) && isfinite(v[2]);
        }

        template< typename Type >
        HOST_DEVICE_PREFIX inline Vector3<Type> abs( const Vector3<Type>& v ) {
            return Vector3<Type>( abs(v[0]), abs(v[1]), abs(v[2]) );
        }

        template< typename Type >
        HOST_DEVICE_PREFIX inline Vector3<Type> fabs( const Vector3<Type>& v ) {
            return Vector3<Type>( fabs(v[0]), fabs(v[1]), fabs(v[2]) );
        }

        template< typename Type >
        HOST_PREFIX std::ostream& operator<<( std::ostream& os, const Vector3<Type>& v )
        {
            return os << "<" << v[0] << "," << v[1] << "," << v[2] << ">";
        }

        template< typename Type >
        HOST_PREFIX std::istream& operator>>( std::istream& is, Vector3<Type>& v )
        {
            if( !is ) return is;

            char bracket1, bracket2, comma1, comma2;
            Type x(0), y(0), z(0);
            const std::istream::pos_type pos( is.tellg() );
            const std::istream::fmtflags oldFlags( is.flags() );

            // Setting the 'skip whitespaces' flag
            is >> std::skipws;

            // Extracting the vector
            if( !(is >> bracket1 >> x >> comma1 >> y >> comma2 >> z >> bracket2) ||
                bracket1 != '<' || comma1 != ',' || comma2 != ',' || bracket2 != '>' ) {
                is.clear();
                is.seekg( pos );
                is.setstate( std::istream::failbit );
                is.flags( oldFlags );
                return is;
            }

            // Transfering the input to the vector values
            v[0] = x; v[1] = y; v[2] = z;

            // Resetting the flags
            is.flags( oldFlags );

            return is;
        }

        //TODO assertion for valid size
        template< typename Buffer_T, typename Type_T >
        HOST_DEVICE_PREFIX Buffer_T * operator<<( Buffer_T * buffer, const Vector3<Type_T>& v ) {

            constexpr auto size = sizeof(Vector3<Type_T>);
            memcpy(buffer, v.data(), sizeof(Vector3<Type_T>));
            return buffer + size;

        }

        //TODO assertion for valid size
        template< typename Buffer_T, typename Type_T >
        HOST_DEVICE_PREFIX Buffer_T * operator>>( Buffer_T * buffer, Vector3<Type_T>& v ) {

            constexpr auto size = sizeof(Vector3<Type_T>);
            memcpy(v.data(), buffer, size);
            return buffer + size;

        }

    } // namespace math

    using math::Vector3;

}

#endif //TURBINECORE_VECTOR3_H
