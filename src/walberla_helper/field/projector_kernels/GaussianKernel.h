
#pragma once

#ifndef TURBINECORE_GAUSSIANKERNEL_H
#define TURBINECORE_GAUSSIANKERNEL_H

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/math/Vector3.h"

namespace turbine_core {

    namespace projectors {

        namespace internal {
            namespace gaussian {
                MANAGED_PREFIX DEVICE_PREFIX Vector3 <real_t> sigma_{};
            }
        }

        template<uint_t widthX, uint_t widthY, uint_t widthZ, bool USE_KCS = false>
        struct GaussianFunction {

        public:

            static constexpr uint_t xWidth{widthX};
            static constexpr uint_t yWidth{widthY};
            static constexpr uint_t zWidth{widthZ};

            HOST_PREFIX static inline void setSigma(const real_t sigma) {
                internal::gaussian::sigma_ = Vector3<real_t>(sigma);
                checkSigma();
            }

            HOST_PREFIX static inline void setSigma(const real_t sigmaX, const real_t sigmaY, const real_t sigmaZ) {
                internal::gaussian::sigma_ = Vector3<real_t>(sigmaX, sigmaY, sigmaZ);
                checkSigma();
            }

            HOST_PREFIX static inline void setSigma(const Vector3 <real_t> &sigma) {
                internal::gaussian::sigma_ = sigma;
                checkSigma();
            }

            HOST_DEVICE_PREFIX static inline real_t getWeight(const real_t &X, const real_t &Y, const real_t &Z,
                                                              const real_t& x, const real_t& y, const real_t& z,
                                                              const real_t &dx = real_t(1), const real_t &dy = real_t(1), const real_t &dz = real_t(1)) {
                return gaussian((X-x) / dx, uint_t(0)) * gaussian((Y-y) / dy, uint_t(1)) * gaussian((Z-z) / dz, uint_t(2));
            }

            HOST_DEVICE_PREFIX static inline real_t getWeight(const Vector3 <real_t> &X, const Vector3 <real_t> &x, const real_t &dx = real_t(1),
                                                              const real_t &dy = real_t(1), const real_t &dz = real_t(1)) {
                return gaussian((X[0]-x[0]) / dx, uint_t(0)) * gaussian((X[1]-x[1]) / dy, uint_t(1)) * gaussian((X[2]-x[2]) / dz, uint_t(2));
            }

        private:

            // TODO hacky solution that might introduce some nasty overhead on GPU due to if condition
            HOST_DEVICE_PREFIX static constexpr uint_t width(const uint_t idx) {
                assert(idx < 3);
                return (idx==0) ? xWidth : ((idx==1) ? yWidth : zWidth);
            }

            template<bool kcs = USE_KCS>
            HOST_DEVICE_PREFIX static inline auto gaussian(const real_t &r, const uint_t &idx)
            -> typename std::enable_if<!kcs, real_t>::type {
                const real_t sigma = internal::gaussian::sigma_[idx];
                assert(sigma > real_t(0.0) && "Sigma must be greater than 0!");
                const real_t rAbs = fabs(r);
                const real_t wHalf = real_t(width(idx)) + real_t(0.5);

                if (rAbs < wHalf) {
                    auto ret = exp(-real_t(0.5) * (rAbs*rAbs / (sigma * sigma))) / (sqrt(real_t(2) * math::pi) * sigma);
                    assert(ret == ret);
                    return ret;
                }
                else
                    return real_t(0);
            }

            template<bool kcs = USE_KCS>
            HOST_DEVICE_PREFIX static inline auto gaussian(const real_t &r, const uint_t &idx)
            -> typename std::enable_if<kcs, real_t>::type {
                const real_t sigma = internal::gaussian::sigma_[idx];
                assert(sigma > real_t(0.0) && "Sigma must be greater than 0!");
                const real_t rAbs = fabs(r);
                const real_t wHalf = real_t(width(idx)) + real_t(0.5);

                if (rAbs < wHalf) {
                    auto ret = exp(real_t(0.5) * (sigma * sigma / (r * r - sigma * sigma) + real_t(1))) / (sqrt(real_t(2) * math::pi) * sigma);
                    assert(ret == ret);
                    return ret;
                }
                else
                    return real_t(0);
            }

            HOST_PREFIX static inline void checkSigma() {
                for(uint_t d = 0; d < 3; ++d) {
                    if (USE_KCS && real_t(width(d)) + real_t(0.5) > internal::gaussian::sigma_[d]) {
                        WALBERLA_ABORT("Kernel width of Gaussian kernel must NOT be greater than 2*sigma for direction " << d << "!")
                    }
                }
            }

        };

        template<bool KCS = false>
        using GaussianFunction3 = GaussianFunction<1,1,1,KCS>;

        template<bool KCS = false>
        using GaussianFunction5 = GaussianFunction<2,2,2,KCS>;

    }

}

#endif //TURBINECORE_GAUSSIANKERNEL_H
