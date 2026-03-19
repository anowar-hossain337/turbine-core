
#pragma once

#ifndef TURBINECORE_MACROSCOPICVARIABLECALCULATOR_H
#define TURBINECORE_MACROSCOPICVARIABLECALCULATOR_H

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/math/Vector3.h"

#include "walberla_helper/Directions.h"


namespace turbine_core {

    namespace field {

        //TODO extend to ForceField = Vector field; currently only multicomponent field possible!

        template< typename Stencil_T, typename PdfField_T, typename ForceField_T, bool zeroCentering = false, bool compressible = false >
        class DensityAndVelocityCalculator {

        public:

            HOST_DEVICE_PREFIX void setForceField( ForceField_T * forceField ) {
                forceField_ = forceField;
            }

            HOST_DEVICE_PREFIX void setPdfField( PdfField_T * pdfField ) {
                pdfField_ = pdfField;
            }

            HOST_DEVICE_PREFIX void get( const Vector3<cell_idx_t> & cell, real_t & density, Vector3<real_t> & velocity ) {
                get(cell[0], cell[1], cell[2], density, velocity);
            }

            HOST_DEVICE_PREFIX void get( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, real_t & density, Vector3<real_t> & velocity ) {

                // reset macroscopic variables
                density = real_t(0);
                velocity = Vector3<real_t>();

                // iterate over stencil
                for( uint_t q = 0; q < Stencil_T::Q; ++q ) {

                    real_t fi = pdfField_->get(x,y,z,q);
                    density += fi;

                    //TODO hack until type conversion in Vector is done
                    Vector3<int> ci_int = stencil::c<Stencil_T>(q);
                    Vector3<real_t> ci{ real_t(ci_int[0]), real_t(ci_int[1]), real_t(ci_int[2]) };
                    velocity += ci * fi;

                }

                // force correction
                if(forceField_ != nullptr) {
                    //TODO assumes dt = 1 -> what about refined/ coarsened regions?
                    Vector3<real_t> force(forceField_->get(x,y,z,0), forceField_->get(x,y,z,1), forceField_->get(x,y,z,2));
                    velocity += real_t(0.5) * force;
                }

                if(zeroCentering) {
                    density += real_t(1);
                }

                if(compressible) {
                    auto rhoInv = real_t(1.0) / density;
                    velocity *= rhoInv;
                }

            }

        private:

            PdfField_T   * pdfField_{nullptr};
            ForceField_T * forceField_{nullptr};

        };

    }

}

#endif //TURBINECORE_MACROSCOPICVARIABLECALCULATOR_H
