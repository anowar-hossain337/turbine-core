
#pragma once

#ifndef TURBINECORE_MATRIX3_H
#define TURBINECORE_MATRIX3_H

#include <iostream>

namespace turbine_core {

    namespace math {

        template< typename T >
        class Matrix3 {

        public:

            HOST_DEVICE_PREFIX Matrix3() {}

            HOST_DEVICE_PREFIX explicit Matrix3( T init ) {
                v_[0] = v_[1] = v_[2] = v_[3] = v_[4] = v_[5] = v_[6] = v_[7] = v_[8] = init;
            }

            HOST_DEVICE_PREFIX Matrix3( const Vector3<T>& a, const Vector3<T>& b, const Vector3<T>& c ) {
                v_[0] = a[0]; v_[1] = b[0]; v_[2] = c[0];
                v_[3] = a[1]; v_[4] = b[1]; v_[5] = c[1];
                v_[6] = a[2]; v_[7] = b[2]; v_[8] = c[2];
            }

            HOST_DEVICE_PREFIX Matrix3( T xx, T xy, T xz,
                                        T yx, T yy, T yz,
                                        T zx, T zy, T zz ) {
                v_[0] = xx; v_[1] = xy; v_[2] = xz;
                v_[3] = yx; v_[4] = yy; v_[5] = yz;
                v_[6] = zx; v_[7] = zy; v_[8] = zz;
            }

            HOST_DEVICE_PREFIX explicit Matrix3( const T* init ) {
                v_[0] = init[0]; v_[1] = init[1]; v_[2] = init[2];
                v_[3] = init[3]; v_[4] = init[4]; v_[5] = init[5];
                v_[6] = init[6]; v_[7] = init[7]; v_[8] = init[8];
            }

            template< typename Axis, typename Angle >
            HOST_DEVICE_PREFIX Matrix3( Vector3<Axis> axis, Angle angle ) {
                const Angle sina( sin(angle) );
                const Angle cosa( cos(angle) );
                const Angle tmp( Angle(1)-cosa );

                normalize(axis);

                v_[0] = cosa + axis[0]*axis[0]*tmp;
                v_[1] = axis[0]*axis[1]*tmp - axis[2]*sina;
                v_[2] = axis[0]*axis[2]*tmp + axis[1]*sina;
                v_[3] = axis[1]*axis[0]*tmp + axis[2]*sina;
                v_[4] = cosa + axis[1]*axis[1]*tmp;
                v_[5] = axis[1]*axis[2]*tmp - axis[0]*sina;
                v_[6] = axis[2]*axis[0]*tmp - axis[1]*sina;
                v_[7] = axis[2]*axis[1]*tmp + axis[0]*sina;
                v_[8] = cosa + axis[2]*axis[2]*tmp;

            }

            HOST_DEVICE_PREFIX Matrix3( const Matrix3& m ) {
                v_[0] = m.v_[0]; v_[1] = m.v_[1]; v_[2] = m.v_[2];
                v_[3] = m.v_[3]; v_[4] = m.v_[4]; v_[5] = m.v_[5];
                v_[6] = m.v_[6]; v_[7] = m.v_[7]; v_[8] = m.v_[8];
            }

            HOST_DEVICE_PREFIX static Matrix3 makeDiagonalMatrix( const T xx, const T yy, const T zz ) {
                return Matrix3( xx, T(), T(),
                                T(),  yy, T(),
                                T(), T(), zz);
            }
            HOST_DEVICE_PREFIX static Matrix3 makeDiagonalMatrix( const T d ) {
                return makeDiagonalMatrix(d,d,d);
            }

            HOST_DEVICE_PREFIX static Matrix3 makeIdentityMatrix() {
                return makeDiagonalMatrix(T(1));
            }

            HOST_DEVICE_PREFIX static Matrix3 makeMatrixFromXYZAngles( const T ax, const T ay, const T az ) {
                return Matrix3<T>(Vector3<T>(1,0,0), ax)
                     * Matrix3<T>(Vector3<T>(0,1,0), ay)
                     * Matrix3<T>(Vector3<T>(0,0,1), az);
            }

            HOST_DEVICE_PREFIX static Matrix3 makeMatrixFromXYZAngles( const Vector3<T> & rotationAngles ) {
                return Matrix3(Vector3<T>(1,0,0), rotationAngles[0])
                     * Matrix3(Vector3<T>(0,1,0), rotationAngles[1])
                     * Matrix3(Vector3<T>(0,0,1), rotationAngles[2]);
            }

            HOST_DEVICE_PREFIX Matrix3& operator= ( T set ) {
                v_[0] = v_[1] = v_[2] = v_[3] = v_[4] = v_[5] = v_[6] = v_[7] = v_[8] = set;
                return *this;
            }

            template< typename Other >
            inline Matrix3& operator=( const Matrix3<Other>& set )
            {
                if(this != &set) {
                    v_[0] = set.v_[0]; v_[1] = set.v_[1]; v_[2] = set.v_[2];
                    v_[3] = set.v_[3]; v_[4] = set.v_[4]; v_[5] = set.v_[5];
                    v_[6] = set.v_[6]; v_[7] = set.v_[7]; v_[8] = set.v_[8];
                }
                return *this;
            }

            HOST_DEVICE_PREFIX Matrix3& operator= ( const Matrix3& set ) {
                if(this != &set) {
                    v_[0] = set.v_[0]; v_[1] = set.v_[1]; v_[2] = set.v_[2];
                    v_[3] = set.v_[3]; v_[4] = set.v_[4]; v_[5] = set.v_[5];
                    v_[6] = set.v_[6]; v_[7] = set.v_[7]; v_[8] = set.v_[8];
                }
                return *this;
            }

            HOST_DEVICE_PREFIX bool operator==( const Matrix3& rhs ) const {
                if( fabs(v_[0]-rhs.v_[0])<1e-8 && fabs(v_[1]-rhs.v_[1])<1e-8 && fabs(v_[2]-rhs.v_[2])<1e-8
                    && fabs(v_[3]-rhs.v_[3])<1e-8 && fabs(v_[4]-rhs.v_[4])<1e-8 && fabs(v_[5]-rhs.v_[5])<1e-8
                    && fabs(v_[6]-rhs.v_[6])<1e-8 && fabs(v_[7]-rhs.v_[7])<1e-8 && fabs(v_[8]-rhs.v_[8])<1e-8)
                    return true;
                return false;
            }

            HOST_DEVICE_PREFIX bool operator!=( const Matrix3& rhs ) const { return !this->operator==(rhs); }

            HOST_DEVICE_PREFIX T& operator[]( uint_t index ) {
                assert(index < 9); return v_[index];
            }
            HOST_DEVICE_PREFIX const T& operator[]( uint_t index ) const {
                assert(index < 9); return v_[index];
            }
            HOST_DEVICE_PREFIX T& operator()( uint_t i, uint_t j ) {
                assert(i < 3 && j < 3); return v_[i*3 + j];
            }
            HOST_DEVICE_PREFIX const T& operator()( uint_t i, uint_t j ) const {
                assert(i < 3 && j < 3); return v_[i*3+j];
            }

            HOST_DEVICE_PREFIX Matrix3& operator+=( const Matrix3& rhs ) {
                v_[0] += rhs.v_[0]; v_[1] += rhs.v_[1]; v_[2] += rhs.v_[2];
                v_[3] += rhs.v_[3]; v_[4] += rhs.v_[4]; v_[5] += rhs.v_[5];
                v_[6] += rhs.v_[6]; v_[7] += rhs.v_[7]; v_[8] += rhs.v_[8];
                return *this;
            }

            HOST_DEVICE_PREFIX Matrix3& operator-=( const Matrix3& rhs ) {
                v_[0] -= rhs.v_[0]; v_[1] -= rhs.v_[1]; v_[2] -= rhs.v_[2];
                v_[3] -= rhs.v_[3]; v_[4] -= rhs.v_[4]; v_[5] -= rhs.v_[5];
                v_[6] -= rhs.v_[6]; v_[7] -= rhs.v_[7]; v_[8] -= rhs.v_[8];
                return *this;
            }

            template<typename Other>
            HOST_DEVICE_PREFIX Matrix3& operator*=( const Matrix3<Other>& rhs ) {
                Matrix3 tmp( v_[0]*rhs.v_[0] + v_[1]*rhs.v_[3] + v_[2]*rhs.v_[6],
                             v_[0]*rhs.v_[1] + v_[1]*rhs.v_[4] + v_[2]*rhs.v_[7],
                             v_[0]*rhs.v_[2] + v_[1]*rhs.v_[5] + v_[2]*rhs.v_[8],
                             v_[3]*rhs.v_[0] + v_[4]*rhs.v_[3] + v_[5]*rhs.v_[6],
                             v_[3]*rhs.v_[1] + v_[4]*rhs.v_[4] + v_[5]*rhs.v_[7],
                             v_[3]*rhs.v_[2] + v_[4]*rhs.v_[5] + v_[5]*rhs.v_[8],
                             v_[6]*rhs.v_[0] + v_[7]*rhs.v_[3] + v_[8]*rhs.v_[6],
                             v_[6]*rhs.v_[1] + v_[7]*rhs.v_[4] + v_[8]*rhs.v_[7],
                             v_[6]*rhs.v_[2] + v_[7]*rhs.v_[5] + v_[8]*rhs.v_[8] );

                return this->operator=( tmp );
            }

            template<typename Other>
            HOST_DEVICE_PREFIX Matrix3 operator+( const Matrix3<Other>& rhs ) const {
                return Matrix3( v_[0] + rhs.v_[0], v_[1] + rhs.v_[1], v_[2] + rhs.v_[2],
                                v_[3] + rhs.v_[3], v_[4] + rhs.v_[4], v_[5] + rhs.v_[5],
                                v_[6] + rhs.v_[6], v_[7] + rhs.v_[7], v_[8] + rhs.v_[8] );
            }

            template<typename Other>
            HOST_DEVICE_PREFIX Matrix3 operator-() const {
                return Matrix3( -v_[0], -v_[1], -v_[2],
                                -v_[3], -v_[4], -v_[5],
                                -v_[6], -v_[7], -v_[8] );
            }

            template<typename Other>
            HOST_DEVICE_PREFIX Matrix3 operator-( const Matrix3<Other>& rhs ) const {
                return Matrix3( v_[0] - rhs.v_[0], v_[1] - rhs.v_[1], v_[2] - rhs.v_[2],
                                v_[3] - rhs.v_[3], v_[4] - rhs.v_[4], v_[5] - rhs.v_[5],
                                v_[6] - rhs.v_[6], v_[7] - rhs.v_[7], v_[8] - rhs.v_[8] );
            }

            HOST_DEVICE_PREFIX Vector3<T> operator* ( const Vector3<T>& rhs ) const {
                return Vector3<T>( v_[0]*rhs[0] + v_[1]*rhs[1] + v_[2]*rhs[2],
                                   v_[3]*rhs[0] + v_[4]*rhs[1] + v_[5]*rhs[2],
                                   v_[6]*rhs[0] + v_[7]*rhs[1] + v_[8]*rhs[2] );
            }

            HOST_DEVICE_PREFIX Matrix3 operator* ( const Matrix3& rhs ) const {
                return Matrix3( v_[0]*rhs.v_[0] + v_[1]*rhs.v_[3] + v_[2]*rhs.v_[6],
                                v_[0]*rhs.v_[1] + v_[1]*rhs.v_[4] + v_[2]*rhs.v_[7],
                                v_[0]*rhs.v_[2] + v_[1]*rhs.v_[5] + v_[2]*rhs.v_[8],
                                v_[3]*rhs.v_[0] + v_[4]*rhs.v_[3] + v_[5]*rhs.v_[6],
                                v_[3]*rhs.v_[1] + v_[4]*rhs.v_[4] + v_[5]*rhs.v_[7],
                                v_[3]*rhs.v_[2] + v_[4]*rhs.v_[5] + v_[5]*rhs.v_[8],
                                v_[6]*rhs.v_[0] + v_[7]*rhs.v_[3] + v_[8]*rhs.v_[6],
                                v_[6]*rhs.v_[1] + v_[7]*rhs.v_[4] + v_[8]*rhs.v_[7],
                                v_[6]*rhs.v_[2] + v_[7]*rhs.v_[5] + v_[8]*rhs.v_[8] );
            }

            template< typename Other >
            HOST_DEVICE_PREFIX Matrix3& operator*=( Other rhs )
            {
                v_[0] *= rhs; v_[1] *= rhs; v_[2] *= rhs;
                v_[3] *= rhs; v_[4] *= rhs; v_[5] *= rhs;
                v_[6] *= rhs; v_[7] *= rhs; v_[8] *= rhs;
                return *this;
            }

            template<typename Other>
            HOST_DEVICE_PREFIX Matrix3 operator*( Other rhs ) const {
                return Matrix3(v_[0] *= rhs, v_[1] *= rhs, v_[2] *= rhs,
                               v_[3] *= rhs, v_[4] *= rhs, v_[5] *= rhs,
                               v_[6] *= rhs, v_[7] *= rhs, v_[8] *= rhs);
            }


            HOST_DEVICE_PREFIX T getDeterminant() const {
                return v_[0]*v_[4]*v_[8] + v_[1]*v_[5]*v_[6] + v_[2]*v_[3]*v_[7] -
                       v_[6]*v_[4]*v_[2] - v_[7]*v_[5]*v_[0] - v_[8]*v_[3]*v_[1];

            }

            HOST_DEVICE_PREFIX Matrix3& transpose() {
                return *this = Matrix3( v_[0], v_[3], v_[6], v_[1], v_[4], v_[7], v_[2], v_[5], v_[8] );
            }

            HOST_DEVICE_PREFIX Matrix3 getTranspose() const {
                return Matrix3( v_[0], v_[3], v_[6], v_[1], v_[4], v_[7], v_[2], v_[5], v_[8] );
            }

            HOST_DEVICE_PREFIX Matrix3& invert() {
                T det = getDeterminant();
                assert( fabs(det-T(0)) > 1e-8 );
                det = T(1) / det;

                return *this = Matrix3( det * ( ( v_[4]*v_[8] ) - ( v_[5]*v_[7] ) ),
                                        det * ( ( v_[7]*v_[2] ) - ( v_[8]*v_[1] ) ),
                                        det * ( ( v_[1]*v_[5] ) - ( v_[2]*v_[4] ) ),
                                        det * ( ( v_[5]*v_[6] ) - ( v_[3]*v_[8] ) ),
                                        det * ( ( v_[8]*v_[0] ) - ( v_[6]*v_[2] ) ),
                                        det * ( ( v_[2]*v_[3] ) - ( v_[0]*v_[5] ) ),
                                        det * ( ( v_[3]*v_[7] ) - ( v_[4]*v_[6] ) ),
                                        det * ( ( v_[6]*v_[1] ) - ( v_[7]*v_[0] ) ),
                                        det * ( ( v_[0]*v_[4] ) - ( v_[1]*v_[3] ) ) );
            }

            HOST_DEVICE_PREFIX Matrix3 getInverse() const {

                T det = getDeterminant();
                assert( fabs(det-T(0)) > 1e-8 );
                det = T(1) / det;

                return Matrix3( det * ( ( v_[4]*v_[8] ) - ( v_[5]*v_[7] ) ),
                                det * ( ( v_[7]*v_[2] ) - ( v_[8]*v_[1] ) ),
                                det * ( ( v_[1]*v_[5] ) - ( v_[2]*v_[4] ) ),
                                det * ( ( v_[5]*v_[6] ) - ( v_[3]*v_[8] ) ),
                                det * ( ( v_[8]*v_[0] ) - ( v_[6]*v_[2] ) ),
                                det * ( ( v_[2]*v_[3] ) - ( v_[0]*v_[5] ) ),
                                det * ( ( v_[3]*v_[7] ) - ( v_[4]*v_[6] ) ),
                                det * ( ( v_[6]*v_[1] ) - ( v_[7]*v_[0] ) ),
                                det * ( ( v_[0]*v_[4] ) - ( v_[1]*v_[3] ) ) );
            }

            HOST_DEVICE_PREFIX bool isSingular() const {
                if( fabs(getDeterminant() - T(0)) < 1e-8 )
                    return true;
                return false;
            }

            HOST_DEVICE_PREFIX bool isSymmetric() const {
                if( fabs(v_[1]-v_[3])<1e-8 && fabs(v_[2]-v_[6])<1e-8 && fabs(v_[5]-v_[7])<1e-8 )
                    return true;
                return false;
            }

            HOST_DEVICE_PREFIX bool isZero() const {
                if(    fabs(v_[0]-T(0))<1e-8 && fabs(v_[1]-T(0))<1e-8 && fabs(v_[2]-T(0))<1e-8
                    && fabs(v_[3]-T(0))<1e-8 && fabs(v_[4]-T(0))<1e-8 && fabs(v_[5]-T(0))<1e-8
                    && fabs(v_[6]-T(0))<1e-8 && fabs(v_[7]-T(0))<1e-8 && fabs(v_[8]-T(0))<1e-8)
                    return true;
                return false;
            }

            HOST_DEVICE_PREFIX T trace() const { return v_[0] + v_[4] + v_[8]; }

            HOST_DEVICE_PREFIX T* data() {return v_;}
            HOST_DEVICE_PREFIX T const * data() const {return v_;}

        private:

            T v_[9]{ T(1), T(0), T(0),
                     T(0), T(1), T(0),
                     T(0), T(0), T(1) };

        };

        template< typename Type >
        HOST_DEVICE_PREFIX bool isnan( const Matrix3<Type>& m ) {
            return isnan(m[0]) && isnan(m[1]) && isnan(m[2]) &&
                   isnan(m[3]) && isnan(m[4]) && isnan(m[5]) &&
                   isnan(m[6]) && isnan(m[7]) && isnan(m[8]);
        }

        template< typename Type >
        HOST_DEVICE_PREFIX bool isinf( const Matrix3<Type>& m ) {
            return isinf(m[0]) && isinf(m[1]) && isinf(m[2]) &&
                   isinf(m[3]) && isinf(m[4]) && isinf(m[5]) &&
                   isinf(m[6]) && isinf(m[7]) && isinf(m[8]);
        }

        template< typename Type >
        HOST_DEVICE_PREFIX bool finite( const Matrix3<Type>& m ) {
            return isfinite(m[0]) && isfinite(m[1]) && isfinite(m[2]) &&
                   isfinite(m[3]) && isfinite(m[4]) && isfinite(m[5]) &&
                   isfinite(m[6]) && isfinite(m[7]) && isfinite(m[8]);
        }

        template< typename Type >
        HOST_DEVICE_PREFIX Matrix3<Type> abs( const Matrix3<Type>& m ) {
            return Matrix3<Type>( abs(m[0]), abs(m[1]), abs(m[2]),
                                  abs(m[3]), abs(m[4]), abs(m[5]),
                                  abs(m[6]), abs(m[7]), abs(m[8]));
        }

        template< typename Type >
        HOST_DEVICE_PREFIX Matrix3<Type> fabs( const Matrix3<Type>& m ) {
            return Matrix3<Type>( fabs(m[0]), fabs(m[1]), fabs(m[2]),
                                  fabs(m[3]), fabs(m[4]), fabs(m[5]),
                                  fabs(m[6]), fabs(m[7]), fabs(m[8]));
        }

        template< typename Type >
        HOST_PREFIX std::ostream& operator<<( std::ostream& os, const Matrix3<Type>& m )
        {
            return os << " ( " << m[0] << " , " << m[1] << " , " << m[2] << " )\n"
                      << " ( " << m[3] << " , " << m[4] << " , " << m[5] << " )\n"
                      << " ( " << m[6] << " , " << m[7] << " , " << m[8] << " )\n";
        }

        //TODO assertion for valid size
        template< typename Buffer_T, typename Type_T >
        HOST_DEVICE_PREFIX Buffer_T * operator<<( Buffer_T * buffer, const Matrix3<Type_T>& m ) {
            constexpr auto size = sizeof(Matrix3<Type_T>);
            memcpy(buffer, m.data(), size);
            return buffer + size;
        }

        //TODO assertion for valid size
        template< typename Buffer_T, typename Type_T >
        HOST_DEVICE_PREFIX Buffer_T * operator>>( Buffer_T * buffer, Matrix3<Type_T>& m ) {
            constexpr auto size = sizeof(Matrix3<Type_T>);
            memcpy(m.data(), buffer, size);
            return buffer + size;
        }

    } // namespace math

    using math::Matrix3;

}

#endif //TURBINECORE_MATRIX3_H
