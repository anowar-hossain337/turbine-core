#pragma once

#ifndef TURBINECORE_DATAINTERPOLATOR_H
#define TURBINECORE_DATAINTERPOLATOR_H

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

#include <cassert>
#include <cmath>

#include <algorithm>

namespace turbine_core {

    namespace data_interpolation {

        template< typename Domain, typename Codomain, bool CUBIC>
        class DataInterpolator {

        public:

            HOST_DEVICE_PREFIX DataInterpolator() {}

            HOST_DEVICE_PREFIX DataInterpolator( const DataInterpolator & interpolator )
                    : nPoints_(interpolator.nPoints_), isWeightSet_(interpolator.isWeightSet_), isSlopeSet_(interpolator.isSlopeSet_)
            {
                x_     = (Domain*)   malloc(nPoints_ * sizeof(Domain)  );
                y_     = (Codomain*) malloc(nPoints_ * sizeof(Codomain));

                for(uint_t i = 0; i < nPoints_; ++i) {
                    x_[i]     = interpolator.x_[i];
                    y_[i]     = interpolator.y_[i];
                }

                if( interpolator.isSlopeSet_ ) {
                    slope_ = (Codomain *) malloc(nPoints_ * sizeof(Codomain));

                    for (uint_t i = 0; i < nPoints_; ++i) {
                        slope_[i] = interpolator.slope_[i];
                    }
                } else if( interpolator.isWeightSet_ ) {
                    weight_ = (Domain *) malloc(nPoints_ * sizeof(Domain));

                    for (uint_t i = 0; i < nPoints_; ++i) {
                        weight_[i] = interpolator.weight_[i];
                    }
                }
            }

            HOST_DEVICE_PREFIX DataInterpolator( DataInterpolator && ) = delete;

            HOST_DEVICE_PREFIX DataInterpolator & operator=( const DataInterpolator & interpolator )
            {
                if(this != &interpolator) {
                    isWeightSet_ = interpolator.isWeightSet_;
                    isSlopeSet_  = interpolator.isSlopeSet_;
                    nPoints_ = interpolator.nPoints_;

                    free(x_); free(y_);
                    x_ = (Domain *) malloc(nPoints_ * sizeof(Domain));
                    y_ = (Codomain *) malloc(nPoints_ * sizeof(Codomain));

                    for (uint_t i = 0; i < nPoints_; ++i) {
                        x_[i] = interpolator.x_[i];
                        y_[i] = interpolator.y_[i];
                    }

                    if ( interpolator.isSlopeSet_ ) {
                        free(slope_);
                        slope_ = (Codomain *) malloc(nPoints_ * sizeof(Codomain));

                        for (uint_t i = 0; i < nPoints_; ++i) {
                            slope_[i] = interpolator.slope_[i];
                        }
                    } else if ( interpolator.isWeightSet_ ) {
                        free(weight_);
                        weight_ = (Domain *) malloc(nPoints_ * sizeof(Domain));

                        for (uint_t i = 0; i < nPoints_; ++i) {
                            weight_[i] = interpolator.weight_[i];
                        }
                    }
                }

                return *this;
            }

            HOST_DEVICE_PREFIX DataInterpolator & operator=( DataInterpolator && ) = delete;

            HOST_DEVICE_PREFIX ~DataInterpolator() {
                free(x_); free(y_); free(weight_); free(slope_);
            }

            HOST_DEVICE_PREFIX uint_t length() const {
                return nPoints_;
            }

            HOST_DEVICE_PREFIX Domain * x() const {
                return x_;
            }

            HOST_DEVICE_PREFIX Codomain * y() const {
                return y_;
            }

            HOST_DEVICE_PREFIX Domain xMin() const {
                return x_[0];
            }

            HOST_DEVICE_PREFIX Domain xMax() const {
                return x_[nPoints_-1];
            }

            /*
             * Based on boost implementation.
             */
            template<bool cubic = CUBIC>
            HOST_DEVICE_PREFIX auto setPoints(const uint_t & nPoints, const Domain * x, const Codomain * y, const uint_t & approximationOrder = 3)
            -> typename std::enable_if<cubic, void>::type {

                using std::min; using std::max;

                assert(nPoints > approximationOrder && "At least 4 data points are required for cubic interpolation.\n");

                nPoints_ = nPoints;

                free(x_); free(y_); free(weight_);
                x_      = (Domain*)   malloc(nPoints_ * sizeof(Domain)  );
                y_      = (Codomain*) malloc(nPoints_ * sizeof(Codomain));
                weight_ = (Codomain*) malloc(nPoints_ * sizeof(Domain));

                for(uint_t i = 0; i < nPoints_; ++i) {
                    x_[i] = x[i];
                    y_[i] = y[i];
                }

                for(uint_t k = 0; k < nPoints_; ++k) {

                    uint_t i_min = max(k - approximationOrder, (uint_t)0);
                    uint_t i_max = k;

                    if( k >= nPoints_ - approximationOrder - 1 ) {
                        i_max = nPoints_ - approximationOrder - 1;
                    }

                    for( uint_t i = i_min; i <= i_max; ++i ) {

                        Domain invProduct{1};
                        uint_t j_max = min( i + approximationOrder, uint_t(nPoints_-1) );

                        for( uint_t j = i; j <= j_max; ++j) {

                            if( j == k )
                                continue;

                            invProduct *= (x_[k] - x_[j]);
                        }

                        if( i % 2 == 0 ) {
                            weight_[k] += Domain(1)/invProduct;
                        } else {
                            weight_[k] -= Domain(1)/invProduct;
                        }
                    }
                }

                isWeightSet_ = true;

            }

            template<bool cubic = CUBIC>
            HOST_DEVICE_PREFIX auto setPoints(const uint_t & nPoints, const Domain * x, const Codomain * y)
            -> typename std::enable_if<!cubic, void>::type {

                assert(nPoints >= 2 && "At least 2 data points are required for linear interpolation.\n");

                nPoints_ = nPoints;

                free(x_); free(y_); free(slope_);
                x_     = (Domain*)   malloc(nPoints_ * sizeof(Domain)  );
                y_     = (Codomain*) malloc(nPoints_ * sizeof(Codomain));
                slope_ = (Codomain*) malloc(nPoints_ * sizeof(Codomain));

                for(uint_t i = 0; i < nPoints_; ++ i) {
                    x_[i] = x[i];
                    y_[i] = y[i];
                }

                for(uint_t i = 0; i < nPoints_ - 1; ++i) {
                    assert(x[i] < x[i+1] && "x data must be strictly monotonously increasing.\n");
                    slope_[i] = (y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i]);
                }

                isSlopeSet_ = true;

                return;
            }

            /*
             * Based on boost implementation.
             */
            template<bool cubic = CUBIC>
            HOST_DEVICE_PREFIX auto operator()(const Domain &x) const
            -> typename std::enable_if<cubic, Codomain>::type {

                Codomain numerator{};
                Domain denominator{};

                for( uint_t i = 0; i < nPoints_; ++i ) {

                    if(x == x_[i]) {
                        return y_[i];
                    }

                    Domain t = weight_[i] / (x - x_[i]);
                    numerator += t * y_[i];
                    denominator += t;

                }

                return numerator / denominator;
            }

            template<bool cubic = CUBIC>
            HOST_DEVICE_PREFIX auto operator()(const Domain &x) const
            -> typename std::enable_if<!cubic, Codomain>::type {

                assert(nPoints_ >= 2 && "At least 2 data points are required for linear interpolation.\n");

                if((x+1e-10) < x_[0]) {
#ifndef NDEBUG
                    printf("In DataInterpolator: x = %f does not lie in data bounds [%f, %f]. Returning lowest y data.\n", x, x_[0], x_[nPoints_-1]);
#endif
                    return y_[0];
                } else if ((x-1e-10) > x_[nPoints_-1]){
#ifndef NDEBUG
                    printf("In DataInterpolator: x = %f does not lie in data bounds [%f, %f]. Returning largest y data.\n", x, x_[0], x_[nPoints_-1]);
#endif
                    return y_[nPoints_-1];
                } else {

                    uint_t lid{0};
                    for(; lid < nPoints_-1; ++lid) {
                        if(x_[lid+1] > x) {
                            break;
                        }
                    }

                    return slope_[lid] * (x - x_[lid]) + y_[lid];

                }

            }

        private:

            uint_t nPoints_{0};

            Domain * x_{nullptr};
            Codomain * y_{nullptr};

            Domain * weight_{nullptr};
            Codomain * slope_{nullptr};

            bool isWeightSet_{false};
            bool isSlopeSet_ {false};

        };


    } // namespace data_interpolation

    namespace math {
        template<typename T>
        class Vector3;
    }

    template<bool CUBIC>
    using ScalarDataInterpolator = data_interpolation::DataInterpolator<real_t, real_t, CUBIC>;

    template<bool CUBIC>
    using Vector3DataInterpolator = data_interpolation::DataInterpolator<real_t, math::Vector3<real_t>, CUBIC>;

} // namespace turbine_core

#endif //TURBINECORE_DATAINTERPOLATOR_H
