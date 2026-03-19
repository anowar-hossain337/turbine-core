
#pragma once

#ifndef TURBINECORE_CONTROLMODEL_RAWSAngles_H
#define TURBINECORE_CONTROLMODEL_RAWSAngles_H

#include "ControlModel.h"
#include "conversion/Conversion.h"

namespace turbine_core {

    namespace control_model {

        class RAWSAnglesModel final: public ControlModel {

            friend class ControlModel;

        public:

            static constexpr Type TYPE = Type::RAWS_ANGLES;

            HOST_DEVICE_PREFIX RAWSAnglesModel() {}

            HOST_PREFIX RAWSAnglesModel(const std::string & controlTable, const Vector3<real_t> initialAngles, const real_t & alphaFilter) {

                alphaFilter_ = alphaFilter;
                previousAngles_ = initialAngles;

                std::ifstream file(controlTable);

                std::string line;

                std::vector<real_t> winds;
                std::vector<Vector3<real_t>> angles;

                if (!file.is_open()) {
                    WALBERLA_ABORT("Error opening file " << controlTable << ", please check file existence.");
                } else {
                    std::getline(file, line); // skip header

                    while (std::getline(file, line)) {

                        bool is_empty = true;
                        for (char ch : line) {
                            is_empty = is_empty && isspace(ch);
                        }
                        if(is_empty)
                            continue;

                        std::stringstream in(line);
                        real_t wind, angleX, angleY, angleZ;
                        in >> wind >> angleX >> angleY >> angleZ;

                        const real_t toRad = walberla::math::pi / real_t(180);

                        winds.push_back(wind * Conversion::C_t() / Conversion::C_l());
                        angles.emplace_back(-angleX * toRad, -angleY * toRad, -angleZ * toRad);
                    }

                    file.close();
                }

                interpolator_.setPoints(winds.size(), winds.data(), angles.data());
            }

            HOST_DEVICE_PREFIX RAWSAnglesModel( const real_t alphaFilter,
                                               const real_t previousRAWS, const real_t angleX,
                                               const real_t angleY, const real_t angleZ, const uint_t nPoints,
                                               real_t * windVelocities, Vector3<real_t> * angles )
                    : alphaFilter_(alphaFilter), previousRAWS_(previousRAWS)
            {
                previousAngles_[0] = angleX;
                previousAngles_[1] = angleY;
                previousAngles_[2] = angleZ;
                interpolator_.setPoints(nPoints, windVelocities, angles);
            }

            HOST_DEVICE_PREFIX RAWSAnglesModel( const RAWSAnglesModel & rawsModel)
                    : interpolator_(rawsModel.interpolator_), previousAngles_(rawsModel.previousAngles_),
                    alphaFilter_(rawsModel.alphaFilter_), previousRAWS_(rawsModel.previousRAWS_) {}

            HOST_DEVICE_PREFIX RAWSAnglesModel & operator=( RAWSAnglesModel && model ) = delete;

            HOST_DEVICE_PREFIX ~RAWSAnglesModel() override {}

            HOST_DEVICE_PREFIX data_interpolation::DataInterpolator<real_t, Vector3<real_t>, false> * interpolator() {
                return &interpolator_;
            }

            HOST_DEVICE_PREFIX auto & alphaFilter() {
                return alphaFilter_;
            }

            HOST_DEVICE_PREFIX auto & previousRAWS() {
                return previousRAWS_;
            }

            HOST_DEVICE_PREFIX auto & previousAngles() {
                return previousAngles_;
            }

        private:

            HOST_DEVICE_PREFIX RAWSAnglesModel * do_clone() const override {
                return new RAWSAnglesModel(*this);
            }

            HOST_DEVICE_PREFIX void do_print() const override {
                printf("RAWSAnglesModel Control Model\n");
            }

            HOST_DEVICE_PREFIX void do_applyControl(Discretisation * discretisation, const uint_t timeStep, 
                                                    const real_t rotorMeanWind, const real_t rotorTorque, 
                                                    const real_t rotorOmega, const Vector3<real_t> rotorWindVaneWind) override {

                previousRAWS_ = alphaFilter_ * previousRAWS_ + real_t(1-alphaFilter_)*rotorMeanWind;

                auto angles = interpolator_(previousRAWS_);
                auto deltaAngles = angles - previousAngles_;
                Quaternion<real_t> quatX(Vector3<real_t>(1,0,0), deltaAngles[0]);
                Quaternion<real_t> quatY(Vector3<real_t>(0,1,0), deltaAngles[1]);
                Quaternion<real_t> quatZ(Vector3<real_t>(0,0,1), deltaAngles[2]);
                auto product = quatX * quatY * quatZ;
                discretisation->updateRelativeOrientation(product);

                previousAngles_ = angles;
            }

            HOST_DEVICE_PREFIX void do_checkNeeds(bool & needsMeanVelocity, bool & needsTorque, bool & needsWindVane) override {
                needsMeanVelocity = true;
                needsTorque = false;
                needsWindVane = false;
                return;
            };

            HOST_DEVICE_PREFIX walberla::mpi::MPISize do_getItemSize() override {
                return 5 * sizeof(real_t);
            };

            HOST_DEVICE_PREFIX void do_fillBuffer(BufferType_T * buffer) override {

                memcpy(buffer, &alphaFilter_, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(buffer, &previousRAWS_, sizeof(real_t));
                buffer = buffer + sizeof(real_t);

                buffer << previousAngles_;
                return;
            };

            HOST_DEVICE_PREFIX void do_unfillBuffer(BufferType_T * buffer) override {

                memcpy(&alphaFilter_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&previousRAWS_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);

                buffer >> previousAngles_;
                return;
            };

            data_interpolation::DataInterpolator<real_t, Vector3<real_t>, false> interpolator_;
            Vector3<real_t> previousAngles_{0, 0., 0.};
            real_t alphaFilter_;
            real_t previousRAWS_{0.};
        };

    } // namespace control_model

} // namespace turbine_core

#endif //TURBINECORE_CONTROLMODEL_RAWSAngles_H
