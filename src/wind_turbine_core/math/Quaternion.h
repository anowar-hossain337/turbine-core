
#pragma once

#ifndef TURBINECORE_QUATERNION_H
#define TURBINECORE_QUATERNION_H

#include <cmath>

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/math/Matrix3.h"

namespace turbine_core {

    namespace math {

        template< typename T >
        class Quaternion {

        public:

            HOST_DEVICE_PREFIX Quaternion() {}

            HOST_DEVICE_PREFIX Quaternion(const T r, const T i, const T j, const T k)
                    : q_{r, i, j, k}
            {
                const T norm = sqrt(r*r + i*i + j*j + k*k);
                q_[0] /= norm;
                q_[1] /= norm;
                q_[2] /= norm;
                q_[3] /= norm;
            }

            template< typename Axis >
            HOST_DEVICE_PREFIX Quaternion(const Vector3 <Axis> & axis, const T angle) {
                const Vector3<Axis> n = normalize(axis);
                const T sina = sin(T(0.5) * angle);

                q_[0] = cos(T(0.5) * angle);
                q_[1] = sina * n[0];
                q_[2] = sina * n[1];
                q_[3] = sina * n[2];
            }

            HOST_DEVICE_PREFIX Quaternion(const T angleX, const T angleY, const T angleZ) {
                reset();
                rotateX( angleX );
                rotateY( angleY );
                rotateZ( angleZ );
            }

            template< typename Other >
            HOST_DEVICE_PREFIX explicit Quaternion(const Vector3 <Other> & euler) {
                reset();
                rotateX( euler[0] );
                rotateY( euler[1] );
                rotateZ( euler[2] );

            }

            HOST_DEVICE_PREFIX Quaternion( const Quaternion& q ) {
                q_[0] = q.q_[0]; q_[1] = q.q_[1];
                q_[2] = q.q_[2]; q_[3] = q.q_[3];
            }

            template<typename Other>
            HOST_DEVICE_PREFIX explicit Quaternion( const Quaternion<Other>& q ) {
                q_[0] = q.q_[0];
                q_[1] = q.q_[1];
                q_[2] = q.q_[2];
                q_[3] = q.q_[3];
            }

            HOST_DEVICE_PREFIX explicit Quaternion( const T * init ) {
                q_[0] = init[0];
                q_[1] = init[1];
                q_[2] = init[2];
                q_[3] = init[3];
            }

            HOST_DEVICE_PREFIX static Quaternion makeQuaternionFromXYZAngles( const T ax, const T ay, const T az) {
                Quaternion<T> q;
                q.reset();
                q.rotateX( ax ); q.rotateY( ay ); q.rotateZ( az );
                return q;
            }

            HOST_DEVICE_PREFIX static Quaternion makeQuaternionFromXYZAngles( const Vector3<T> & angles ) {
                return makeQuaternionFromXYZAngles(angles[0], angles[1], angles[2]);
            }

            HOST_DEVICE_PREFIX Quaternion& operator= ( const Quaternion& rhs ) {
                if(this != &rhs) {
                    q_[0] = rhs.q_[0];
                    q_[1] = rhs.q_[1];
                    q_[2] = rhs.q_[2];
                    q_[3] = rhs.q_[3];
                }
                return *this;
            }

            template< typename Other >
            HOST_DEVICE_PREFIX Quaternion& operator= ( const Quaternion<Other>& rhs ) {
                if(this != &rhs) {
                    q_[0] = rhs.q_[0];
                    q_[1] = rhs.q_[1];
                    q_[2] = rhs.q_[2];
                    q_[3] = rhs.q_[3];
                }
                return *this;
            }

            HOST_DEVICE_PREFIX T& operator[]( size_t index ) {
                assert( index < 4 );
                return q_[index];
            }

            HOST_DEVICE_PREFIX T operator[]( size_t index ) const {
                assert( index < 4 );
                return q_[index];
            }

            HOST_DEVICE_PREFIX Quaternion& set( const T r, const T i, const T j, const T k ) {
                assert( fabs( r*r + i*i + j*j + k*k - T(1) ) < 1e-8 );
                q_[0] = r;
                q_[1] = i;
                q_[2] = j;
                q_[3] = k;
                return *this;
            }

            HOST_DEVICE_PREFIX void reset() {
                q_[0] = T(1);
                q_[1] = q_[2] = q_[3] = T(0);
            }

            HOST_DEVICE_PREFIX Quaternion operator*(const Quaternion & q) const {
                const T w = q_[0] * q.q_[0] - q_[1] * q.q_[1] - q_[2] * q.q_[2] - q_[3] * q.q_[3];
                const T x = q_[0] * q.q_[1] + q.q_[0] * q_[1] + q_[2] * q.q_[3] - q_[3] * q.q_[2];
                const T y = q_[0] * q.q_[2] + q.q_[0] * q_[2] + q_[3] * q.q_[1] - q_[1] * q.q_[3];
                const T z = q_[0] * q.q_[3] + q.q_[0] * q_[3] + q_[1] * q.q_[2] - q_[2] * q.q_[1];

                return Quaternion{w,x,y,z};
            }

            HOST_DEVICE_PREFIX Quaternion& invert() {
                q_[1] = -q_[1];
                q_[2] = -q_[2];
                q_[3] = -q_[3];
                return *this;
            }

            HOST_DEVICE_PREFIX Quaternion getInverse() const {
                return Quaternion( q_[0], -q_[1], -q_[2], -q_[3] );
            }

            HOST_DEVICE_PREFIX Matrix3<T> toRotationMatrix() const {
                const T xx = T(1) - T(2) * ( q_[2] * q_[2] + q_[3] * q_[3] );
                const T yy = T(1) - T(2) * ( q_[1] * q_[1] + q_[3] * q_[3] );
                const T zz = T(1) - T(2) * ( q_[1] * q_[1] + q_[2] * q_[2] );

                const T xy = T(2) * ( q_[1] * q_[2] - q_[0] * q_[3] );
                const T yx = T(2) * ( q_[1] * q_[2] + q_[0] * q_[3] );
                const T xz = T(2) * ( q_[1] * q_[3] + q_[0] * q_[2] );
                const T zx = T(2) * ( q_[1] * q_[3] - q_[0] * q_[2] );
                const T yz = T(2) * ( q_[2] * q_[3] - q_[0] * q_[1] );
                const T zy = T(2) * ( q_[2] * q_[3] + q_[0] * q_[1] );

                return Matrix3<T>( xx, xy, xz,
                                   yx, yy, yz,
                                   zx, zy, zz);
            }

            HOST_DEVICE_PREFIX T getAngle() const {
                return T(2) * atan2(sqrt(T(1) - q_[0] * q_[0]), q_[0]);
            }

            HOST_DEVICE_PREFIX Vector3<T> getAxis() const {
                const T normsq = T(1) - q_[0] * q_[0];
                if( fabs(normsq) < 1e-8 )
                    return Vector3<T>(1,0,0);
                const T norm = sqrt(normsq);
                return Vector3<T>(q_[1], q_[2], q_[3]) / norm;
            }

            HOST_DEVICE_PREFIX void rotateX( T angle ) {
                const T sina( sin( angle*T(0.5) ) );
                const T cosa( cos( angle*T(0.5) ) );

                const Quaternion q( cosa*q_[0] - sina*q_[1],
                                    cosa*q_[1] + sina*q_[0],
                                    cosa*q_[2] - sina*q_[3],
                                    cosa*q_[3] + sina*q_[2] );

                this->operator=( q );
            }

            HOST_DEVICE_PREFIX void rotateY( T angle ) {
                const T sina( sin( angle*T(0.5) ) );
                const T cosa( cos( angle*T(0.5) ) );

                const Quaternion q( cosa*q_[0] - sina*q_[2],
                                    cosa*q_[1] + sina*q_[3],
                                    cosa*q_[2] + sina*q_[0],
                                    cosa*q_[3] - sina*q_[1] );

                this->operator=( q );
            }

            HOST_DEVICE_PREFIX void rotateZ( T angle ) {
                const T sina( sin( angle*T(0.5) ) );
                const T cosa( cos( angle*T(0.5) ) );

                const Quaternion q( cosa*q_[0] - sina*q_[3],
                                    cosa*q_[1] - sina*q_[2],
                                    cosa*q_[2] + sina*q_[1],
                                    cosa*q_[3] + sina*q_[0] );

                this->operator=( q );
            }

            HOST_DEVICE_PREFIX void swap( Quaternion& q ) {
                Quaternion tmp = q;
                q = *this;
                *this = tmp;
            }

            HOST_DEVICE_PREFIX Vector3<T> getEulerAnglesXYZ() const {
                Vector3<T> eulerAngles;

                eulerAngles[0] = atan2( T(2) * ( q_[0] * q_[1] + q_[2] * q_[3] ), T(1) - T(2) * ( q_[1] * q_[1] + q_[2] * q_[2] ) );
                eulerAngles[1] = asin ( T(2) * ( q_[0] * q_[2] - q_[3] * q_[1] ) );
                eulerAngles[2] = atan2( T(2) * ( q_[0] * q_[3] + q_[1] * q_[2] ), T(1) - T(2) * ( q_[2] * q_[2] + q_[3] * q_[3] ) );

                return eulerAngles;
            }


            HOST_DEVICE_PREFIX T * data() {return q_;}
            HOST_DEVICE_PREFIX T const * data() const {return q_;}

            /*template< typename Other >
            HOST_DEVICE_PREFIX Vector3<T> rotate( const Vector3<Other>& v ) const {
                const auto w( q_[1]*v[0] + q_[2]*v[1] + q_[3]*v[2] );
                const auto x( q_[0]*v[0] - q_[3]*v[1] + q_[2]*v[2] );
                const auto y( q_[0]*v[1] - q_[1]*v[2] + q_[3]*v[0] );
                const auto z( q_[0]*v[2] - q_[2]*v[0] + q_[1]*v[1] );

                return Vector3<T>( q_[0]*x + q_[1]*w + q_[2]*z - q_[3]*y,
                                   q_[0]*y + q_[2]*w + q_[3]*x - q_[1]*z,
                                   q_[0]*z + q_[3]*w + q_[1]*y - q_[2]*x );

            }*/

            template< typename Other >
            HOST_DEVICE_PREFIX Vector3<T> rotate( const Vector3<Other>& v ) const {
                const Vector3<T> r {q_[1], q_[2], q_[3]};
                return v + T(2.0) * r % ( r % v + q_[0] * v );
            }


        private:

            T q_[4]{ T(1), T(0), T(0), T(0) };

        };

        //TODO assertion for valid size
        template< typename Buffer_T, typename Type_T >
        HOST_DEVICE_PREFIX Buffer_T * operator<<( Buffer_T * buffer, const Quaternion<Type_T>& q ) {
            constexpr auto size = sizeof(Quaternion<Type_T>);
            memcpy(buffer, q.data(), size);
            return buffer + size;
        }

        //TODO assertion for valid size
        template< typename Buffer_T, typename Type_T >
        HOST_DEVICE_PREFIX Buffer_T * operator>>( Buffer_T * buffer, Quaternion<Type_T>& q ) {
            constexpr auto size = sizeof(Quaternion<Type_T>);
            memcpy(q.data(), buffer, size);
            return buffer + size;
        }

    }

    using math::Quaternion;

}

#endif //TURBINECORE_QUATERNION_H
