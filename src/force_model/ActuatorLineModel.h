#pragma once

#ifndef TURBINECORE_ACTUATORLINEMODEL_H
#define TURBINECORE_ACTUATORLINEMODEL_H

#include "ForceModel.h"
#include "TipLossModel.h"

#include <fstream>
#include <vector>

#include <core/math/Constants.h>

#include "conversion/Conversion.h"

#include "wind_turbine_core/math/Matrix2.h"
#include "wind_turbine_core/math/Vector2.h"

namespace turbine_core {

    namespace force_model {

        template< typename Interpolator_T, typename ForceDistributor_T, typename TipLossModel_T >
        class ActuatorLineModel final: public ForceModel {

        public:

            HOST_DEVICE_PREFIX ActuatorLineModel( const real_t elementLength, const uint_t nPoints, const ActuatorData * data, TipLossModel_T tiplossModel )
                    : nPoints_(nPoints), elementLength_(elementLength), tiplossModel_(tiplossModel)
            {
                points_ = (ActuatorData*) malloc(nPoints_ * sizeof(ActuatorData));
                processed_ = (bool*) malloc(2 * nPoints_ * sizeof(bool));

                for (uint_t i = 0; i < nPoints_; ++i) {
                    new (static_cast<ActuatorData*>(points_) + i) ActuatorData(data[i]);
                    processed_[2 * i] = false;
                    processed_[2 * i + 1] = false;
                }

                macroProcessed_ = processed_;
                remotePoint_ = processed_ + nPoints_;
            }

            template<bool CUBIC>
            HOST_DEVICE_PREFIX ActuatorLineModel( const uint_t nPoints, const uint_t nInterpolationPoints, const real_t * span, const aerodynamics::AirfoilPolar * airfoilPolars, const ScalarDataInterpolator<CUBIC> & widthInterpolator, TipLossModel_T tiplossModel )
                    : nPoints_(nPoints)
            {
                points_ = (ActuatorData*) malloc(nPoints_ * sizeof(ActuatorData));
                processed_ = (bool*) malloc(2 * nPoints_ * sizeof(bool));
                elementLength_ = (span[nInterpolationPoints-1] - span[0]) / real_t(nPoints);

                for (unsigned int i = 0; i < nPoints; ++i) {

                    const real_t tmpSpan = real_t(i + 0.5) * elementLength_;

                    uint_t lid{0};
                    for(; lid < nInterpolationPoints-1; ++lid) {
                        if(tmpSpan < span[lid+1])
                            break;
                    }

                    real_t weight = (tmpSpan - span[lid]) / (span[lid+1] - span[lid]);

                    aerodynamics::AirfoilPolar polar(weight, airfoilPolars[lid], airfoilPolars[lid+1]);

                    new (static_cast<ActuatorData*>(points_) + i) ActuatorData(widthInterpolator(tmpSpan), polar);

                    processed_[2 * i] = false;
                    processed_[2 * i + 1] = false;

                }

                macroProcessed_ = processed_;
                remotePoint_ = processed_ + nPoints_;

                tiplossModel_ = tiplossModel;
            }

            HOST_DEVICE_PREFIX ActuatorLineModel(const ActuatorLineModel & model)
                    : nPoints_(model.nPoints_), elementLength_(model.elementLength_), tiplossModel_(model.tiplossModel_)
            {
                points_ = (ActuatorData*) malloc(nPoints_ * sizeof(ActuatorData));
                processed_ = (bool*) malloc(2 * nPoints_ * sizeof(bool));

                for (uint_t i = 0; i < nPoints_; ++i) {
                    new (static_cast<ActuatorData*>(points_) + i) ActuatorData(model.points_[i]);
                    processed_[2 * i] = model.processed_[2 * i];
                    processed_[2 * i + 1] = model.processed_[2 * i + 1];
                }

                macroProcessed_ = processed_;
                remotePoint_ = processed_ + nPoints_;
            }

            HOST_DEVICE_PREFIX ActuatorLineModel( ActuatorLineModel && model ) = delete;

            HOST_DEVICE_PREFIX ActuatorLineModel & operator=(const ActuatorLineModel & model) {

                if(this != &model) {

                    nPoints_ = model.nPoints_;

                    free(points_);
                    points_ = (ActuatorData*) malloc(nPoints_ * sizeof(ActuatorData));

                    free(processed_);
                    processed_ = (bool*) malloc(2 * nPoints_ * sizeof(bool));

                    for (uint_t i = 0; i < nPoints_; ++i) {
                        new (static_cast<ActuatorData*>(points_) + i) ActuatorData(model.points_[i]);
                        processed_[2 * i] = model.processed_[2 * i];
                        processed_[2 * i + 1] = model.processed_[2 * i + 1];
                    }

                    macroProcessed_ = processed_;
                    remotePoint_ = processed_ + nPoints_;

                    tiplossModel_ = model.tiplossModel_;
                }

                return *this;
            }

            HOST_DEVICE_PREFIX ActuatorLineModel & operator=( ActuatorLineModel && model ) = delete;

            HOST_DEVICE_PREFIX ~ActuatorLineModel() override {
                for(uint_t i = 0; i < nPoints_; ++i) {
                    (static_cast<ActuatorData*>(points_) + i)->~ActuatorData();
                }
                free(points_);
                free(processed_);
            }

            HOST_DEVICE_PREFIX auto & elementLength() {
                return elementLength_;
            }

            HOST_DEVICE_PREFIX ActuatorData * points() override {
                return points_;
            }

            HOST_DEVICE_PREFIX bool * remotePoint() override {
                return remotePoint_;
            }

        private:

            HOST_DEVICE_PREFIX ActuatorLineModel * do_clone() const override {
                return new ActuatorLineModel(*this);
            }

            HOST_DEVICE_PREFIX void do_print() const override {
                printf("ActuatorLine Model nPoints = %lu\n", nPoints_);

                for (uint_t i = 0; i < nPoints_; ++i) {
//                    printf("ActuatorLine Model data[%lu].density = %f\n", i, points_[i].density);
                }

            }

            HOST_DEVICE_PREFIX void do_evaluateDensityAndVelocity( const blockforest::BlockInfo & blockInfo,
                                                                   const projectors::Projectors & projectors,
                                                                   const Point3 <real_t> *physicalPoints,
                                                                   uint_t pointIdy) override {

                /**
                 * Convention for coordinate systems:
                 * - physicalPoints is expressed completely in global coordinates
                 *      -> point.velocity is global
                 *      -> point.orientation transforms from local to global coordinates
                 **/

                auto densityInterpolator  = projectors.getDensityInterpolator<Interpolator_T>();
                auto velocityInterpolator  = projectors.getVelocityInterpolator<Interpolator_T>();

                const auto & blockBB = blockInfo.getAABB();

                if (pointIdy < nPoints_) {
                    const auto &point = physicalPoints[pointIdy];
                    const auto &currentPosition = point.position;

                    if (blockBB.contains(currentPosition)) {

                        macroProcessed_[pointIdy] = true;

                        real_t localDensity{};
                        real_t velocity[3]{};
                        densityInterpolator->get(currentPosition[0], currentPosition[1], currentPosition[2],
                                                 &localDensity);
                        velocityInterpolator->get(currentPosition[0], currentPosition[1], currentPosition[2], velocity);
                        Vector3 <real_t> windVelocityGlobal{velocity[0], velocity[1], velocity[2]};

                        Vector3 <real_t> prevWindVelocityGlobal = points_[pointIdy].windVelocityGlobal;

                        const auto particleVel = windVelocityGlobal.length();
                        const auto smoothing = real_t(exp(-2.0 * math::pi * particleVel / points_[pointIdy].width * 1.0));
                        windVelocityGlobal = smoothing * prevWindVelocityGlobal + (real_t(1.0) - smoothing) * windVelocityGlobal;

                        points_[pointIdy].windVelocityGlobal = windVelocityGlobal;

                        assert(localDensity > real_t(0.0) && "interpolated density must be larger than 0!");

                        points_[pointIdy].density = localDensity;
                        auto relativeVelocityLocal = point.orientation.getInverse().rotate(windVelocityGlobal - point.velocity);

                        assert(!isnan(relativeVelocityLocal) && "component of interpolated velocity is nan!");

                        points_[pointIdy].velocityLocal = relativeVelocityLocal;

                    }
                }// loop nPoints

            }

            HOST_DEVICE_PREFIX void do_calculateForces( const Point3<real_t> * physicalPoints,
                                                        const Point3<real_t> * parentPoint ) override {

                for( uint_t i = 0; i < nPoints_; ++i ) {

                    //TODO also calculate remote points here? -> get rid of force communication?
                    // point not calculated, i.e. not valid!
                    if( !(macroProcessed_[i] || remotePoint_[i]) ) {
                        points_[i].angleOfAttack = real_t(0);
                        points_[i].forcesGlobal = Vector3<real_t>{};
                        continue;
                    }

                    const auto & point = physicalPoints[i];

                    auto relativeVelocityLocal = points_[i].velocityLocal;

                    // neglect span-wise velocity in rel. velocity calc. and calculate velocity magnitude
                    relativeVelocityLocal[2] = real_t(0.0);
                    real_t velocityMagnitude = relativeVelocityLocal.length();

                    auto aoa = real_t(atan2( relativeVelocityLocal[0], relativeVelocityLocal[1] ));
                    points_[i].angleOfAttack = aoa;

                    auto polars = points_[i].polar(aoa);

                    auto coefficients = Matrix2<real_t>{  real_t(cos(aoa)),  real_t(sin(aoa)),
                                                          -real_t(sin(aoa)), +real_t(cos(aoa)) } * polars;

                    // blade element local x-axis is pointing upward -> minus sign
                    const real_t tipLoss = tiplossModel_.evaluateTipLoss(parentPoint, point, physicalPoints[0], physicalPoints[nPoints_-1], elementLength_, relativeVelocityLocal);
                    assert(!(tipLoss != tipLoss) && "calculated tiploss is nan!");
                    auto forcesLocal = coefficients * real_t(0.5) * tipLoss * points_[i].density * real_t(pow(velocityMagnitude, 2.0))
                                  * points_[i].width * elementLength_;

                    assert(!isnan(forcesLocal) && "component of calculated force is nan!");

                    Vector3<real_t> forcesLocal3d{};
                    forcesLocal3d[0] = forcesLocal[0];
                    forcesLocal3d[1] = forcesLocal[1];
                    forcesLocal3d[2] = real_t(0);

                    // rotate forces into global reference frame
                    points_[i].forcesGlobal = point.orientation.rotate(forcesLocal3d);

                } // loop nPoints

            }

            HOST_DEVICE_PREFIX void do_spreadForces( const blockforest::BlockInfo & blockInfo,
                                                     const projectors::Projectors & projectors,
                                                     const Point3 <real_t> *physicalPoints,
                                                     uint_t pointIdy) override {

                auto forceDistributor = projectors.getForceDistributor<ForceDistributor_T>();

                using Kernel_T = typename ForceDistributor_T::Kernel;
                const auto gl = projectors.nGhostLayers();
                // points that distribute into ghost layers (might be outside of ghost layers)
                const Vector3<real_t> kernelWidth{
                    gl + Kernel_T::xWidth,
                    gl + Kernel_T::yWidth,
                    gl + Kernel_T::zWidth
                };

                const auto & blockBB = blockInfo.getAABB().getExtended(kernelWidth);

                if (pointIdy < nPoints_) {

                    const auto &point = physicalPoints[pointIdy];

                    if ((macroProcessed_[pointIdy] || remotePoint_[pointIdy]) && blockBB.contains(point.position)) {
                        forceDistributor->distribute(point.position[0], point.position[1], point.position[2],
                                                     (-points_[pointIdy].forcesGlobal).data());
                    }

                } // loop nPoints

            }

            HOST_DEVICE_PREFIX void do_resetCommunicationData() override {
                memset(processed_, 0, 2 * nPoints_ * sizeof(bool));
            }

        public:

            uint_t nPoints_{0};
            ActuatorData * points_{nullptr};

            bool * processed_{nullptr};
            bool * macroProcessed_{nullptr};
            bool * remotePoint_{nullptr};

            real_t elementLength_{};

            TipLossModel_T tiplossModel_{};
        };

    }

}

#endif //TURBINECORE_ACTUATORLINEMODEL_H
