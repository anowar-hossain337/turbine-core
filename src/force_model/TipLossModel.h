#pragma once

#ifndef TURBINECORE_TIPLOSSMODEL_H
#define TURBINECORE_TIPLOSSMODEL_H

#include <core/math/Constants.h>

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "point3/Point3.h"

namespace turbine_core {

    namespace force_model {

        class TipLossModel {

        public:

            HOST_DEVICE_PREFIX TipLossModel() {}

            HOST_DEVICE_PREFIX virtual ~TipLossModel() {}

            HOST_DEVICE_PREFIX TipLossModel * clone() const {
                return do_clone();
            }

            HOST_DEVICE_PREFIX real_t evaluateTipLoss(const Point3<real_t> * parentPoint,
                                                const Point3<real_t> & point,
                                                const Point3<real_t> & startPoint,
                                                const Point3<real_t> & endPoint,
                                                real_t& pointLength,
                                                const Vector3<real_t> & relativeVelocity) {
                return do_evaluateTipLoss(parentPoint, point, startPoint, endPoint, pointLength, relativeVelocity);
            }

        private:

            HOST_DEVICE_PREFIX virtual TipLossModel * do_clone() const = 0;

            HOST_DEVICE_PREFIX virtual real_t do_evaluateTipLoss(const Point3<real_t> * parentPoint,
                                                         const Point3<real_t> & point,
                                                         const Point3<real_t>  & startPoint,
                                                         const Point3<real_t>  & endPoint,
                                                         real_t& pointLength,
                                                         const Vector3<real_t> & relativeVelocity) = 0;
        };


        class NoTipLoss final: public TipLossModel {

            friend class TipLossModel;

        public:

            HOST_DEVICE_PREFIX NoTipLoss() {}

            HOST_DEVICE_PREFIX NoTipLoss( const NoTipLoss & ) {}

            HOST_DEVICE_PREFIX NoTipLoss( NoTipLoss && model ) = delete;

            HOST_DEVICE_PREFIX NoTipLoss & operator=(const NoTipLoss &) {
                return *this;
            }

            HOST_DEVICE_PREFIX NoTipLoss & operator=( NoTipLoss && model ) = delete;

            HOST_DEVICE_PREFIX ~NoTipLoss() override {}

        private:

            HOST_DEVICE_PREFIX NoTipLoss * do_clone() const override {
                return new NoTipLoss(*this);
            }

            HOST_DEVICE_PREFIX real_t do_evaluateTipLoss(const Point3<real_t> * parentPoint,
                                                 const Point3<real_t> & point,
                                                 const Point3<real_t> & startPoint,
                                                 const Point3<real_t> & endPoint,
                                                 real_t& pointLength,
                                                 const Vector3<real_t> & relativeVelocity) override {
                return 1.;
            }
        };


        class PrandtlTipLoss final: public TipLossModel {

            friend class TipLossModel;

        public:

            HOST_DEVICE_PREFIX PrandtlTipLoss() {}

            HOST_DEVICE_PREFIX PrandtlTipLoss( const PrandtlTipLoss & ) {}

            HOST_DEVICE_PREFIX PrandtlTipLoss( NoTipLoss && model ) = delete;

            HOST_DEVICE_PREFIX PrandtlTipLoss & operator=(const PrandtlTipLoss &) {
                return *this;
            }

            HOST_DEVICE_PREFIX PrandtlTipLoss & operator=( PrandtlTipLoss && model ) = delete;

            HOST_DEVICE_PREFIX ~PrandtlTipLoss() override {}

        private:
            HOST_DEVICE_PREFIX PrandtlTipLoss * do_clone() const override {
                return new PrandtlTipLoss(*this);
            }

            HOST_DEVICE_PREFIX real_t do_evaluateTipLoss(const Point3<real_t> * parentPoint,
                                                 const Point3<real_t> & point,
                                                 const Point3<real_t>  & startPoint,
                                                 const Point3<real_t>  & endPoint,
                                                 real_t& pointLength,
                                                 const Vector3<real_t> & relativeVelocity) override {

                auto radius = (point.position - parentPoint->position).length();

                auto relVelHub = parentPoint->orientation.getInverse().rotate(point.orientation.rotate(relativeVelocity));
                auto inPlaneVelocity = sqrt(relVelHub[0]*relVelHub[0] + relVelHub[1]*relVelHub[1]);
                auto phi = -atan2( relVelHub[2], inPlaneVelocity );

                auto hubRadius = (startPoint.position-parentPoint->position).length() - real_t(.5)*pointLength;
                auto rotorRadius = (endPoint.position-parentPoint->position).length() + real_t(.5)*pointLength;

                //TODO: use real blade number instead of 3
                //TODO: implement hubLoss
                auto tipLoss = real_t(2.)/math::pi*acos(exp(-real_t(3.)*(rotorRadius-radius)/(real_t(2.)*radius*sin(phi))));
                return tipLoss;
            }
        };
    }

    using force_model::TipLossModel;

}

#endif //TURBINECORE_TIPLOSSMODEL_H
