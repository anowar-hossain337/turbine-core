
#pragma once

#ifndef TURBINECORE_ACTUATORDISKMODEL_H
#define TURBINECORE_ACTUATORDISKMODEL_H

#include <core/math/Constants.h>
#include <core/math/all.h>

#include "wind_turbine_core/ProjectDefines.h"

#include "force_model/ForceModel.h"

namespace turbine_core {

    namespace force_model {

        template< typename Interpolator_T, typename ForceDistributor_T >
        class ActuatorDiskModel final: public ForceModel {

        public:

            HOST_DEVICE_PREFIX explicit ActuatorDiskModel( const real_t radius, const real_t drag )
                    : radius_{fmax(radius, real_t(1))}, drag_{drag}
            {
                point_ = new ActuatorData;
            }

            HOST_DEVICE_PREFIX ActuatorDiskModel(const ActuatorDiskModel & model)
                    : radius_(model.radius_), drag_(model.drag_)
            {
                if(model.point_)
                    point_ = new ActuatorData(*model.point_);
            }

            HOST_DEVICE_PREFIX ActuatorDiskModel( ActuatorDiskModel && model ) = delete;

            HOST_DEVICE_PREFIX ActuatorDiskModel & operator=(const ActuatorDiskModel & model) {
                if(this != &model) {
                    radius_ = model.radius_;
                    drag_ = model.drag_;
                    point_ = new ActuatorData(*model.point_);
                }
                return *this;
            }

            HOST_DEVICE_PREFIX ActuatorDiskModel & operator=( ActuatorDiskModel && model ) = delete;

            HOST_DEVICE_PREFIX ~ActuatorDiskModel() override {
                delete point_;
            }

            HOST_DEVICE_PREFIX ActuatorData * points() override {
                return point_;
            }

            HOST_DEVICE_PREFIX bool * remotePoint() override {
                return nullptr;
            }

        private:

            HOST_DEVICE_PREFIX ActuatorDiskModel * do_clone() const override {
                return new ActuatorDiskModel(*this);
            }

            HOST_DEVICE_PREFIX void do_print() const override {

                printf("ActuatorDiskModel radius and drag = %f %f\n", radius_, drag_);

            }

            HOST_DEVICE_PREFIX void do_evaluateDensityAndVelocity( const blockforest::BlockInfo & blockInfo,
                                                                   const projectors::Projectors & projectors,
                                                                   const Point3 <real_t> *physicalPoints,
                                                                   uint_t pointIdy) override {

                if (pointIdy > 0) return;
                const auto & domainBB = blockInfo.getAABB();

                const auto & domainBBMin = domainBB.min();
                assert( domainBBMin[0] >= 0 && domainBBMin[1] >= 0 && domainBBMin[2] >= 0 && "Min corner of domain must not have negative coordinates." );

                auto & point = physicalPoints[0];

                auto currentPosition = point.position;
                if( domainBB.contains(currentPosition) ) {

                    auto densityInterpolator  = projectors.getDensityInterpolator<Interpolator_T>();
                    auto velocityInterpolator  = projectors.getVelocityInterpolator<Interpolator_T>();

                    real_t localDensity{};
                    real_t velocity[3]{};
                    densityInterpolator->get(currentPosition[0], currentPosition[1], currentPosition[2], &localDensity);
                    velocityInterpolator->get(currentPosition[0], currentPosition[1], currentPosition[2], velocity);
                    Vector3<real_t> windVelocityGlobal{velocity[0], velocity[1], velocity[2]};

                    point_->density = localDensity;

                    const auto refLength = uint_t(ceil(radius_));

                    Vector3<real_t> prevWindVelocityGlobal = point_->windVelocityGlobal;
                    const auto particleVel = prevWindVelocityGlobal.length();

                    const auto smoothing = math::min(real_t(exp( - 2.0 * math::pi * particleVel / refLength * 1.0 )), 0.98);
                    point_->windVelocityGlobal = smoothing * prevWindVelocityGlobal + (real_t(1.0) - smoothing) * windVelocityGlobal;
                }
            }

            HOST_DEVICE_PREFIX void do_calculateForces( const Point3<real_t> * physicalPoints,
                                                        const Point3<real_t> * ) override {

                Vector3<real_t> windVelocityLocal = physicalPoints[0].orientation.getInverse().rotate(point_->windVelocityGlobal);
                real_t velocityMagnitude = fabs(windVelocityLocal[2]);
                auto fz = real_t(drag_ * real_t(0.5) * point_->density * pow(velocityMagnitude, 2.0));

                Vector3<real_t> forcesLocal{0, 0, fz};
                point_->forcesGlobal = physicalPoints[0].orientation.rotate(forcesLocal);
            }

            HOST_DEVICE_PREFIX void do_spreadForces( const blockforest::BlockInfo & blockInfo,
                                                     const projectors::Projectors & projectors,
                                                     const Point3 <real_t> *physicalPoints,
                                                     uint_t pointIdy) override {

                auto forceDistributor = projectors.getForceDistributor<ForceDistributor_T>();

                using K_T = typename ForceDistributor_T::Kernel;
                const Vector3<real_t> kernelWidth (K_T::xWidth, K_T::yWidth, K_T::zWidth);

                const auto & extendedBlockBB = blockInfo.getAABB().getExtended(kernelWidth);

                const uint_t diameter = 2u*uint_t(ceil(radius_));

                const auto & diskCenter = physicalPoints[0].position;
                const auto & orientation = physicalPoints[0].orientation;

                auto min = Vector3<real_t>(1, 1, 1) * INFINITY;
                auto max = Vector3<real_t>(-1, -1, -1) * INFINITY;

                // Corners of aligned disk bounding box
                Vector3<real_t> corners[8];
                corners[0] = Vector3<real_t>{real_t(-0.5), -radius_, -radius_};
                corners[1] = Vector3<real_t>{real_t(-0.5), -radius_, +radius_};
                corners[2] = Vector3<real_t>{real_t(-0.5), +radius_, -radius_};
                corners[3] = Vector3<real_t>{real_t(-0.5), +radius_, +radius_};
                corners[4] = Vector3<real_t>{real_t(+0.5), -radius_, -radius_};
                corners[5] = Vector3<real_t>{real_t(+0.5), -radius_, +radius_};
                corners[6] = Vector3<real_t>{real_t(+0.5), +radius_, -radius_};
                corners[7] = Vector3<real_t>{real_t(+0.5), +radius_, +radius_};

                for (uint_t i = 0; i < 8; ++i) {
                    auto transformed = orientation.getInverse().rotate(corners[i]);
                    for (uint_t d = 0; d < 3; ++d) {
                        min[d] = ::fmin(min[d], transformed[d]);
                        max[d] = ::fmax(max[d], transformed[d]);
                    }
                }

                //FIXME the explicit coordinate exchange (x-z) accounts for the difference in local coordinate system, i.e.,
                // the global x-axis is the disk-local z-axis; can this be automised so that it works for any disks ?
                const int xMin{floor(min[2])}, xMax{ceil(max[2])}, yMin{floor(min[1])}, yMax{ceil(max[1])}, zMin{floor(min[0])}, zMax{ceil(max[0])};

//                printf("xMin = %d, \txMax = %d, \nyMin = %d, \tyMax = %d \nzMin = %d, \tzMax = %d\n\n", xMin, xMax, yMin, yMax, zMin, zMax);

//                for( int k = zMin; k < zMax; ++k ) {
//                    for( int j = yMin; j < yMax; ++j ) {
//                        for( int i = xMin; i < xMax; ++i ) {
                int dimX = xMax - xMin;
                int dimY = yMax - yMin;
                int dimZ = zMax - zMin;

                if (pointIdy >= dimX * dimY * dimZ)
                    return;

                int k = pointIdy / (dimX * dimY) + zMin;
                int j = (pointIdy % (dimX * dimY)) / dimX + yMin;
                int i = (pointIdy % dimX) + xMin;

                if (i < xMax && j < yMax && k < zMax) {
                    // TODO: improve this temporary fix: an epsilon is introduced to avoid a
                    // "-0" distance after projection into the disk point reference frame, which
                    // excludes this valid point from the if-condition.
                    Vector3 <real_t> currentPosition{
                            diskCenter[0] + real_t(i) + real_t(1e-5),
                            diskCenter[1] + real_t(j),
                            diskCenter[2] + real_t(k)
                    };

                    if (!extendedBlockBB.contains(currentPosition)) {
                        return;
                    }

                    auto &point = physicalPoints[0];
                    Vector3 <real_t> distance = currentPosition - point.position;
                    auto distanceLocal = point.orientation.getInverse().rotate(distance);

                    if (distanceLocal[0] * distanceLocal[0] + distanceLocal[1] * distanceLocal[1] <= radius_ * radius_
                        && distanceLocal[2] < 1. && distanceLocal[2] >= 0.) {
                        forceDistributor->distribute(currentPosition[0], currentPosition[1], currentPosition[2],
                                                     point_->forcesGlobal.data());
                    }
                }
            }

        public:

            real_t radius_{};
            real_t drag_{0};
            ActuatorData * point_{nullptr};
        };

    }

}

#endif //TURBINECORE_ACTUATORDISKMODEL_H
