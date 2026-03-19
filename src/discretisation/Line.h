
#pragma once

#ifndef TURBINECORE_LINE_H
#define TURBINECORE_LINE_H

#include "Discretisation.h"

#include <core/DataTypes.h>

#include "wind_turbine_core/ProjectDefines.h"

#include "data_interpolator/DataInterpolator.h"

namespace turbine_core {

    namespace discretisation {

        class Line final : public Discretisation {

            typedef Point3<real_t> PointType;

        public:

            HOST_DEVICE_PREFIX Line() {}

            HOST_DEVICE_PREFIX explicit Line(const uint_t & nPoints )
                    : nPoints_(nPoints)
            {
                free(points_);
                points_ = (PointType*)malloc(nPoints_ * sizeof(PointType));
            }

            HOST_DEVICE_PREFIX Line(const Vector3<real_t> & start, const Vector3<real_t> & end )
                    : start_{}, end_{}, nPoints_{0}, points_{nullptr}
            {
                start_.position = start;
                end_.position = end;
            }

            HOST_DEVICE_PREFIX Line(const PointType & start, const PointType & end,
                                    const uint_t & nPoints, PointType * points)
                    : start_(start), end_(end), nPoints_(nPoints)
            {
                free(points_);
                points_ = (PointType*)malloc(nPoints_ * sizeof(PointType));

                for( uint_t i = 0; i < nPoints_; ++i ) {
                    points_[i] = points[i];
                }
            }

            HOST_DEVICE_PREFIX Line( const Vector3<real_t> & relativePosition, const Vector3<real_t> & relativeTranslationalVelocity,
                                     const Vector3<real_t> & relativeRotationalVelocity, const Quaternion<real_t> & relativeOrientation )
                    : relativePosition_(relativePosition), relativeTranslationVelocity_(relativeTranslationalVelocity),
                      relativeRotationalVelocity_(relativeRotationalVelocity), relativeOrientation_(relativeOrientation),
                      rotationMatrix_(Quaternion<real_t>::makeQuaternionFromXYZAngles(relativeRotationalVelocity_)),
                      fullRotationMatrix_(rotationMatrix_.getInverse()), //NOTE to account for initial update, where no rotation should be performed
                      startRotationMatrix_(rotationMatrix_.getInverse()),
                      totalTranslation_(-relativeTranslationVelocity_) //TODO transform in correct cooridnate system!
            {}

            HOST_DEVICE_PREFIX Line(const Line & line)
                    : start_(line.start_), end_(line.end_), nPoints_(line.nPoints_), elementLength_(line.elementLength_),
                      relativePosition_(line.relativePosition_), relativeTranslationVelocity_(line.relativeTranslationVelocity_),
                      relativeRotationalVelocity_(line.relativeRotationalVelocity_), relativeOrientation_(line.relativeOrientation_),
                      rotationMatrix_(line.rotationMatrix_), fullRotationMatrix_(line.fullRotationMatrix_), startRotationMatrix_(line.startRotationMatrix_),
                      totalTranslation_(line.totalTranslation_)
            {
                free(points_);
                if(nPoints_ != 0) {
                    points_ = (PointType*)malloc(line.nPoints_ * sizeof(PointType));

                            for( uint_t i = 0; i < nPoints_; ++i ) {
                              points_[i] = line.points_[i];
                            }
                } else {
                    points_ = nullptr;
                }

            }

            HOST_DEVICE_PREFIX Line(const Line & line, const uint_t & nPoints, PointType * points) {

                if(this != &line) {

                    start_ = line.start_;
                    end_ = line.end_;

                    nPoints_ = nPoints;
                    elementLength_ = line.elementLength_;

                    free(points_);
                    if(nPoints != 0) {
                                points_ = (PointType*)malloc(nPoints * sizeof(PointType));

                    for( uint_t i = 0; i < nPoints; ++i ) {
                        points_[i] = points[i];
                                }
                    } else {
                    points_ = nullptr;
                    }

                    relativePosition_ = line.relativePosition_;
                    relativeTranslationVelocity_ = line.relativeTranslationVelocity_;
                    relativeRotationalVelocity_ = line.relativeRotationalVelocity_;
                    relativeOrientation_ = line.relativeOrientation_;
                    rotationMatrix_ = line.rotationMatrix_;
                    fullRotationMatrix_ = line.fullRotationMatrix_;
                    startRotationMatrix_ = line.startRotationMatrix_;
                    totalTranslation_ = line.totalTranslation_;

                }

            }

            HOST_DEVICE_PREFIX Line(Line && line ) = delete;

            HOST_DEVICE_PREFIX Line & operator=(const Line & line) {

                if(this != &line) {

                    start_ = line.start_;
                    end_ = line.end_;

                    nPoints_ = line.nPoints_;
                    elementLength_ = line.elementLength_;

                    free(points_);

                    if(line.nPoints_ != 0) {
                        points_ = (PointType *)malloc(line.nPoints_ * sizeof(PointType));

                        for (uint_t i = 0; i < nPoints_; ++i) {
                            points_[i] = line.points_[i];
                        }
                    } else {
                        points_ = nullptr;
                    }

                    relativePosition_ = line.relativePosition_;
                    relativeTranslationVelocity_ = line.relativeTranslationVelocity_;
                    relativeRotationalVelocity_ = line.relativeRotationalVelocity_;
                    relativeOrientation_ = line.relativeOrientation_;
                    rotationMatrix_ = line.rotationMatrix_;
                    totalTranslation_ = line.totalTranslation_;
                    startRotationMatrix_ = line.startRotationMatrix_;

                }

                return *this;
            }

            HOST_DEVICE_PREFIX Line & operator=(Line && Line ) = delete;

            HOST_DEVICE_PREFIX ~Line() override {
                free(points_);
            }

            template<bool CUBIC>
            HOST_DEVICE_PREFIX void setPoints( const uint_t nPoints,
                                               const Vector3DataInterpolator<CUBIC> * positionInterpolator,
                                               const Vector3DataInterpolator<CUBIC> * angleInterpolator ) {

                nPoints_ = nPoints;
                free(points_);
                points_ = (PointType*)malloc(nPoints_ * sizeof(PointType));

                const uint_t nInterpolationPoints = positionInterpolator->length();

                start_.relativePosition = positionInterpolator->y()[0];
                end_.relativePosition = positionInterpolator->y()[nInterpolationPoints-1];

                // convention: pitch etc around negative local z-axis
                start_.relativeOrientation = Quaternion<real_t>::makeQuaternionFromXYZAngles(-angleInterpolator->y()[0]);
                end_.relativeOrientation = Quaternion<real_t>::makeQuaternionFromXYZAngles(-angleInterpolator->y()[nInterpolationPoints-1]);

                auto span = positionInterpolator->x();
                elementLength_ = (span[nInterpolationPoints-1] - span[0]) / real_t(nPoints);

                for(uint_t i = 0; i < nPoints; ++i) {

                    const real_t tmpSpan = real_t(i + 0.5) * elementLength_;
                    points_[i] = PointType();

                    points_[i].relativePosition = positionInterpolator->operator()(tmpSpan);
                    points_[i].relativeOrientation = Quaternion<real_t>::makeQuaternionFromXYZAngles(-angleInterpolator->operator()(tmpSpan));

                }

            }

            HOST_DEVICE_PREFIX void calculateAzimuthAngle( const Point3<real_t> * reference, const uint_t idx, real_t & azimuthAngle ) const {

                //FIXME: the current implementation does not work for a VAWT, a more generic approach is requiered

                // Point to nacelle distance in nacelle reference frame
                Vector3 <real_t> pointDistance = reference->orientation.getInverse().rotate(points_[idx].position - reference->position);

                // Azimuth angle evaluation
                azimuthAngle =
                        std::atan2(pointDistance[0], pointDistance[1]) / math::pi *
                        real_t(180);

                // 0-360 bounds
                azimuthAngle = fmod(azimuthAngle + 360., real_t(360));
            }

            HOST_DEVICE_PREFIX PointType * points() override {
                return points_;
            }

            HOST_DEVICE_PREFIX uint_t nPoints() override {
                return nPoints_;
            }

            HOST_DEVICE_PREFIX auto startPoint() {
                return start_;
            }

            HOST_DEVICE_PREFIX auto endPoint() {
                return end_;
            }

        private:

            HOST_DEVICE_PREFIX Line * do_clone() const override {
                return new Line(*this);
            }

            HOST_DEVICE_PREFIX void do_print() const override {

                printf("PolyLine relativePosition =    [%f, %f, %f]\n", relativePosition_[0], relativePosition_[1], relativePosition_[2]);

                printf("PolyLine start   = [%f, %f, %f]\n", start_.position[0], start_.position
                [1], start_.position[2]);

                for(uint_t i = 0; i < nPoints_; ++i) {
                    printf("PolyLine point %lu = [%f, %f, %f]\n", i, points_[i].position[0], points_[i].position[1], points_[i].position[2]);
                }
                printf("PolyLine   end   = [%f, %f, %f]\n", end_.position[0],   end_.position[1],   end_.position[2]);
            }

            HOST_DEVICE_PREFIX void do_update( const Point3<real_t> * reference, const uint_t ) override {

                fullRotationMatrix_ = fullRotationMatrix_ * rotationMatrix_;

                auto globalRelativePosition = reference->orientation.rotate(relativePosition_);
                auto globalRelativeTranslationalVelocity = start_.orientation.rotate(relativeTranslationVelocity_);

                auto rotInc = relativeRotationalVelocity_;
                 startRotationMatrix_ = startRotationMatrix_ * Quaternion<real_t>::makeQuaternionFromXYZAngles(rotInc);

                start_.orientation = reference->orientation * relativeOrientation_ * startRotationMatrix_;
                start_.position    = reference->position + globalRelativePosition + reference->orientation.rotate(start_.relativePosition);

                // update velocities //TODO what to do in first step?
                totalTranslation_ = totalTranslation_ + globalRelativeTranslationalVelocity;

                start_.rotationalVelocity = reference->rotationalVelocity + start_.orientation.rotate(relativeRotationalVelocity_);

                start_.velocity = reference->velocity + globalRelativeTranslationalVelocity
                                  + reference->rotationalVelocity % globalRelativePosition;

                start_.position = start_.position + totalTranslation_;

                // update points depending on start
                end_.updatePoint(start_, fullRotationMatrix_, relativeTranslationVelocity_);

                for(uint_t i = 0; i < nPoints_; ++i) {
                    points_[i].updatePoint(start_, fullRotationMatrix_);
                }

            }

            HOST_DEVICE_PREFIX Point3<real_t> * do_getReferencePointForChild() override {
                return &end_;
            }

            HOST_DEVICE_PREFIX void do_setRelativeRotationalVelocity( const Vector3<real_t> & relativeRotationalVelocity) override {
                relativeRotationalVelocity_ = relativeRotationalVelocity;
            }

            HOST_DEVICE_PREFIX void do_setRelativeTranslationVelocity( const Vector3<real_t> & relativeTranslationVelocity) override {
                relativeTranslationVelocity_ = relativeTranslationVelocity;
            }

            HOST_DEVICE_PREFIX void do_updateRelativeOrientation( const Quaternion<real_t> & relativeOrientation) override {
                relativeOrientation_ = relativeOrientation_ * relativeOrientation;
            }

            HOST_DEVICE_PREFIX const Quaternion<real_t> & do_getRelativeOrientation() const override {
                return rotationMatrix_;
            }

            uint_t nPoints_{0};
            PointType * points_{nullptr};

        public:

            PointType start_{};
            PointType end_{};

            real_t elementLength_{};

            Vector3<real_t> relativePosition_{};
            Vector3<real_t> relativeTranslationVelocity_{};
            Vector3<real_t> relativeRotationalVelocity_{};

            Quaternion<real_t> relativeOrientation_{};

            Quaternion<real_t> rotationMatrix_{};
            Quaternion<real_t> fullRotationMatrix_{};
            Quaternion<real_t> startRotationMatrix_{};

            Vector3<real_t> totalTranslation_{};

        };
    }

}

#endif //TURBINECORE_LINE_H
