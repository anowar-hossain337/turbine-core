
#pragma once

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/ProjectDefines.h"

#include <cassert>

#include <core/debug/Debug.h>
#include <core/Abort.h>
#include <core/logging/Logging.h>

namespace turbine_core {

    namespace conversion {

        struct ConversionFactors {

        public:

            HOST_PREFIX static void calculateConversionFactors(real_t length_SI, real_t length_LU, real_t velocity_SI,
                                                               real_t velocity_LU, real_t density_SI);

            HOST_PREFIX static real_t C_l();
            HOST_PREFIX static real_t C_t();
            HOST_PREFIX static real_t C_m();

            HOST_PREFIX static void print();

        private:

            static bool calculated_;
            static real_t factors_[3];

        };

        bool ConversionFactors::calculated_ = false;
        real_t ConversionFactors::factors_[3] = {real_t(-1), real_t(-1), real_t(-1)};

        HOST_PREFIX void ConversionFactors::calculateConversionFactors(const real_t length_SI, const real_t length_LU,
                                                                       const real_t velocity_SI, const real_t velocity_LU,
                                                                       const real_t density_SI) {

            factors_[0] = length_SI / length_LU;
            factors_[1] = velocity_LU / velocity_SI * factors_[0];
            factors_[2] = density_SI / 1. * factors_[0] * factors_[0] * factors_[0];

            WALBERLA_CHECK_GREATER(factors_[0], real_t(0.0),
                                    "Conversion factor calculation failed for: C_l");
            WALBERLA_CHECK_GREATER(factors_[1], real_t(0.0),
                                    "Conversion factor calculation failed for: C_t");
            WALBERLA_CHECK_GREATER(factors_[2], real_t(0.0),
                                    "Conversion factor calculation failed for: C_m");

            calculated_ = true;

        }

        HOST_PREFIX real_t ConversionFactors::C_l() {

            assert(calculated_ && "Please calculate the conversion factors before accessing!");
            return factors_[0];
        }

        HOST_PREFIX real_t ConversionFactors::C_t() {
            assert(calculated_ && "Please calculate the conversion factors before accessing!");
            return factors_[1];
        }

        HOST_PREFIX real_t ConversionFactors::C_m() {
            assert(calculated_ && "Please calculate the conversion factors before accessing!");
            return factors_[2];
        }

        HOST_PREFIX void ConversionFactors::print() {

            WALBERLA_LOG_INFO_ON_ROOT("C_l = " << C_l() << "\n" <<
                                      "C_t = " << C_t() << "\n" <<
                                      "C_m = " << C_m())

        }

    } // namespace conversion

    using Conversion = conversion::ConversionFactors;

} // namespace turbine_core