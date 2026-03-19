
#pragma once

#ifndef TURBINECORE_MATRIX2_H
#define TURBINECORE_MATRIX2_H

#include <iostream>

namespace turbine_core {

    namespace math {

        template< typename T >
        class Matrix2 {

        public:

            HOST_DEVICE_PREFIX Matrix2() {}

            HOST_DEVICE_PREFIX explicit Matrix2( T init ) {
                v_[0] = v_[1] = v_[2] = v_[3] = init;
            }

            HOST_DEVICE_PREFIX Matrix2( const Vector2<T>& a, const Vector2<T>& b ) {
                v_[0] = a[0]; v_[1] = b[0];
                v_[2] = a[0]; v_[3] = b[1];
            }

            HOST_DEVICE_PREFIX Matrix2( T xx, T xy,
                                        T yx, T yy ) {
                v_[0] = xx; v_[1] = xy;
                v_[2] = yx; v_[3] = yy;
            }

            HOST_DEVICE_PREFIX explicit Matrix2( const T* init ) {
                v_[0] = init[0]; v_[1] = init[1];
                v_[2] = init[2]; v_[3] = init[3];
            }

            HOST_DEVICE_PREFIX Matrix2( const Matrix2& m ) {
                v_[0] = m.v_[0]; v_[1] = m.v_[1];
                v_[2] = m.v_[2]; v_[3] = m.v_[3];
            }

            HOST_DEVICE_PREFIX static Matrix2 makeDiagonalMatrix( const T xx, const T yy ) {
                return Matrix2(  xx, T(),
                                T(), yy );
            }
            HOST_DEVICE_PREFIX static Matrix2 makeDiagonalMatrix( const T d ) {
                return makeDiagonalMatrix(d,d);
            }

            HOST_DEVICE_PREFIX static Matrix2 makeIdentityMatrix() {
                return makeDiagonalMatrix(T(1));
            }

            HOST_DEVICE_PREFIX static Matrix2 makeMatrixFromXYZAngles( const T ax, const T ay ) {
                return Matrix2<T>(Vector2<T>(1,0), ax)
                     * Matrix2<T>(Vector2<T>(0,1), ay);
            }

            HOST_DEVICE_PREFIX static Matrix2 makeMatrixFromXYZAngles( const Vector2<T> & rotationAngles ) {
                return Matrix2(Vector2<T>(1,0), rotationAngles[0])
                     * Matrix2(Vector2<T>(0,1), rotationAngles[1]);
            }

            HOST_DEVICE_PREFIX static Matrix2 makeRotationMatrix( T angle ) {
                const T sina( sin(angle) );
                const T cosa( cos(angle) );

                return Matrix2<T>(cosa, -sina,
                                  sina,  cosa);
            }

            HOST_DEVICE_PREFIX Matrix2& operator= ( T set ) {
                v_[0] = v_[1] = v_[2] = v_[3] = v_[4] = v_[5] = v_[6] = v_[7] = v_[8] = set;
                return *this;
            }

            template< typename Other >
            inline Matrix2& operator=( const Matrix2<Other>& set ) {
                if(this != &set) {
                    v_[0] = set.v_[0]; v_[1] = set.v_[1];
                    v_[2] = set.v_[2]; v_[3] = set.v_[3];
                }
                return *this;
            }

            HOST_DEVICE_PREFIX Matrix2& operator= ( const Matrix2& set ) {
                if(this != &set) {
                    v_[0] = set.v_[0]; v_[1] = set.v_[1];
                    v_[2] = set.v_[2]; v_[3] = set.v_[3];
                }
                return *this;
            }

            HOST_DEVICE_PREFIX bool operator==( const Matrix2& rhs ) const {
                if( fabs(v_[0]-rhs.v_[0])<1e-8 && fabs(v_[1]-rhs.v_[1])<1e-8 &&
                    fabs(v_[2]-rhs.v_[2])<1e-8 && fabs(v_[3]-rhs.v_[3])<1e-8   )
                    return true;
                return false;
            }

            HOST_DEVICE_PREFIX bool operator!=( const Matrix2& rhs ) const { return !this->operator==(rhs); }

            HOST_DEVICE_PREFIX T& operator[]( uint_t index ) {
                assert(index < 4); return v_[index];
            }
            HOST_DEVICE_PREFIX const T& operator[]( uint_t index ) const {
                assert(index < 4); return v_[index];
            }
            HOST_DEVICE_PREFIX T& operator()( uint_t i, uint_t j ) {
                assert(i < 2 && j < 2); return v_[i*2 + j];
            }
            HOST_DEVICE_PREFIX const T& operator()( uint_t i, uint_t j ) const {
                assert(i < 2 && j < 2); return v_[i*2+j];
            }

            HOST_DEVICE_PREFIX Matrix2& operator+=( const Matrix2& rhs ) {
                v_[0] += rhs.v_[0]; v_[1] += rhs.v_[1];
                v_[2] += rhs.v_[2]; v_[3] += rhs.v_[3];
                return *this;
            }

            HOST_DEVICE_PREFIX Matrix2& operator-=( const Matrix2& rhs ) {
                v_[0] -= rhs.v_[0]; v_[1] -= rhs.v_[1];
                v_[2] -= rhs.v_[2]; v_[3] -= rhs.v_[3];
                return *this;
            }

            template<typename Other>
            HOST_DEVICE_PREFIX Matrix2& operator*=( const Matrix2<Other>& rhs ) {
                Matrix2 tmp( v_[0]*rhs.v_[0] + v_[1]*rhs.v_[2],
                             v_[0]*rhs.v_[1] + v_[1]*rhs.v_[3],
                             v_[2]*rhs.v_[0] + v_[3]*rhs.v_[2],
                             v_[2]*rhs.v_[1] + v_[3]*rhs.v_[3] );

                return this->operator=( tmp );
            }

            template<typename Other>
            HOST_DEVICE_PREFIX Matrix2 operator+( const Matrix2<Other>& rhs ) const {
                return Matrix2( v_[0] + rhs.v_[0], v_[1] + rhs.v_[1],
                                v_[2] + rhs.v_[2], v_[3] + rhs.v_[3] );
            }

            template<typename Other>
            HOST_DEVICE_PREFIX Matrix2 operator-() const {
                return Matrix2( -v_[0], -v_[1],
                                -v_[2], -v_[3] );
            }

            template<typename Other>
            HOST_DEVICE_PREFIX Matrix2 operator-( const Matrix3<Other>& rhs ) const {
                return Matrix2( v_[0] - rhs.v_[0], v_[1] - rhs.v_[1],
                                v_[2] - rhs.v_[2], v_[3] - rhs.v_[3] );
            }

            HOST_DEVICE_PREFIX Vector2<T> operator* ( const Vector2<T>& rhs ) const {
                return Vector2<T>( v_[0]*rhs[0] + v_[1]*rhs[1],
                                   v_[2]*rhs[0] + v_[3]*rhs[1] );
            }

            HOST_DEVICE_PREFIX Matrix2 operator* ( const Matrix2& rhs ) const {
                return Matrix2( v_[0]*rhs.v_[0] + v_[1]*rhs.v_[2],
                                v_[0]*rhs.v_[1] + v_[1]*rhs.v_[3],
                                v_[2]*rhs.v_[0] + v_[3]*rhs.v_[2],
                                v_[2]*rhs.v_[1] + v_[3]*rhs.v_[3] );
            }

            template< typename Other >
            HOST_DEVICE_PREFIX Matrix2& operator*=( Other rhs )
            {
                v_[0] *= rhs; v_[1] *= rhs;
                v_[2] *= rhs; v_[3] *= rhs;
                return *this;
            }

            template<typename Other>
            HOST_DEVICE_PREFIX Matrix2 operator*( Other rhs ) const {
                return Matrix2(v_[0] *= rhs, v_[1] *= rhs,
                               v_[2] *= rhs, v_[3] *= rhs );
            }


            HOST_DEVICE_PREFIX T getDeterminant() const {
                return v_[0]*v_[3] - v_[1]*v_[3];

            }

            HOST_DEVICE_PREFIX Matrix2& transpose() {
                return *this = Matrix2( v_[0], v_[2], v_[1], v_[3] );
            }

            HOST_DEVICE_PREFIX Matrix2 getTranspose() const {
                return Matrix2( v_[0], v_[2], v_[1], v_[3] );
            }

            HOST_DEVICE_PREFIX Matrix2& invert() {
                T det = getDeterminant();
                assert( fabs(det-T(0)) > 1e-8 );
                det = T(1) / det;

                return *this = Matrix2(   det * v_[3], - det * v_[1],
                                        - det * v_[2],   det * v_[0] );
            }

            HOST_DEVICE_PREFIX Matrix2 getInverse() const {

                T det = getDeterminant();
                assert( fabs(det-T(0)) > 1e-8 );
                det = T(1) / det;

                return Matrix2(   det * v_[3], - det * v_[1],
                                - det * v_[2],   det * v_[0] );
            }

            HOST_DEVICE_PREFIX bool isSingular() const {
                if( fabs(getDeterminant() - T(0)) < 1e-8 )
                    return true;
                return false;
            }

            HOST_DEVICE_PREFIX bool isSymmetric() const {
                if( fabs(v_[1]-v_[2])<1e-8 )
                    return true;
                return false;
            }

            HOST_DEVICE_PREFIX bool isZero() const {
                if(    fabs(v_[0]-T(0))<1e-8 && fabs(v_[1]-T(0))<1e-8 &&
                       fabs(v_[2]-T(0))<1e-8 && fabs(v_[3]-T(0))<1e-8   )
                    return true;
                return false;
            }

            HOST_DEVICE_PREFIX T trace() const { return v_[0] + v_[3]; }

            HOST_DEVICE_PREFIX T* data() {return v_;}
            HOST_DEVICE_PREFIX T const * data() const {return v_;}

        private:

            T v_[4]{ T(1), T(0),
                     T(0), T(1) };

        };

        template< typename Type >
        HOST_DEVICE_PREFIX bool isnan( const Matrix2<Type>& m ) {
            return isnan(m[0]) && isnan(m[1]) &&
                   isnan(m[2]) && isnan(m[3]);
        }

        template< typename Type >
        HOST_DEVICE_PREFIX bool isinf( const Matrix2<Type>& m ) {
            return isinf(m[0]) && isinf(m[1]) &&
                   isinf(m[2]) && isinf(m[3]);
        }

        template< typename Type >
        HOST_DEVICE_PREFIX bool finite( const Matrix2<Type>& m ) {
            return isfinite(m[0]) && isfinite(m[1]) &&
                   isfinite(m[2]) && isfinite(m[3]);
        }

        template< typename Type >
        HOST_DEVICE_PREFIX Matrix2<Type> abs( const Matrix2<Type>& m ) {
            return Matrix2<Type>( abs(m[0]), abs(m[1]),
                                  abs(m[2]), abs(m[3]) );
        }

        template< typename Type >
        HOST_DEVICE_PREFIX Matrix2<Type> fabs( const Matrix2<Type>& m ) {
            return Matrix2<Type>( fabs(m[0]), fabs(m[1]),
                                  fabs(m[2]), fabs(m[3]));
        }

        template< typename Type >
        HOST_PREFIX std::ostream& operator<<( std::ostream& os, const Matrix2<Type>& m )
        {
            return os << " ( " << m[0] << " , " << m[1]  << " )\n"
                      << " ( " << m[2] << " , " << m[3]  << " )\n";
        }

        //TODO assertion for valid size
        template< typename Buffer_T, typename Type_T >
        HOST_DEVICE_PREFIX Buffer_T * operator<<( Buffer_T * buffer, const Matrix2<Type_T>& m ) {
            constexpr auto size = sizeof(Matrix2<Type_T>);
            memcpy(buffer, m.data(), size);
            return buffer + size;
        }

        //TODO assertion for valid size
        template< typename Buffer_T, typename Type_T >
        HOST_DEVICE_PREFIX Buffer_T * operator>>( Buffer_T * buffer, Matrix2<Type_T>& m ) {
            constexpr auto size = sizeof(Matrix2<Type_T>);
            memcpy(m.data(), buffer, size);
            return buffer + size;
        }

} // namespace math

using math::Matrix2;

}

#endif //TURBINECORE_MATRIX2_H
