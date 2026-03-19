
#pragma once

#ifndef TURBINECORE_CONTROLMODEL_RAWSOmega_H
#define TURBINECORE_CONTROLMODEL_RAWSOmega_H

#include "ControlModel.h"
#include "conversion/Conversion.h"

namespace turbine_core {

    namespace control_model {

        class RAWSOmegaModel final: public ControlModel {

            friend class ControlModel;

        public:

            static constexpr Type TYPE = Type::RAWS_OMEGA;

            HOST_DEVICE_PREFIX RAWSOmegaModel() {}

            HOST_PREFIX RAWSOmegaModel(const std::string & controlTable, const real_t & alphaFilter) {

                alphaFilter_ = alphaFilter;

                std::ifstream file(controlTable);

                std::string line;

                std::vector<real_t> winds;
                std::vector<Vector3<real_t>> omegas;

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
                        real_t wind, omegaX, omegaY, omegaZ;
                        in >> wind >> omegaX >> omegaY >> omegaZ;

                        winds.push_back(wind * Conversion::C_t() / Conversion::C_l());
                        omegas.emplace_back(-omegaX * Conversion::C_t(), -omegaY * Conversion::C_t(), -omegaZ * Conversion::C_t());
                    }

                    file.close();
                }

                interpolator_.setPoints(winds.size(), winds.data(), omegas.data());
            }

            HOST_DEVICE_PREFIX RAWSOmegaModel( const real_t alphaFilter,
                                               const real_t previousRAWS, const uint_t nPoints,
                                               real_t * windVelocities, Vector3<real_t> * rotationalVelocities )
                    : alphaFilter_(alphaFilter), previousRAWS_(previousRAWS)
            {
                interpolator_.setPoints(nPoints, windVelocities, rotationalVelocities);
            }

            HOST_DEVICE_PREFIX RAWSOmegaModel( const RAWSOmegaModel & rawsModel)
                : interpolator_(rawsModel.interpolator_), alphaFilter_(rawsModel.alphaFilter_), previousRAWS_(rawsModel.previousRAWS_) {}

            HOST_DEVICE_PREFIX RAWSOmegaModel & operator=( RAWSOmegaModel && model ) = delete;

            HOST_DEVICE_PREFIX ~RAWSOmegaModel() override = default;

            HOST_DEVICE_PREFIX data_interpolation::DataInterpolator<real_t, Vector3<real_t>, false> * interpolator() {
                return &interpolator_;
            }

            HOST_DEVICE_PREFIX auto & alphaFilter() {
                return alphaFilter_;
            }

            HOST_DEVICE_PREFIX auto & previousRAWS() {
                return previousRAWS_;
            }

        private:

            HOST_DEVICE_PREFIX RAWSOmegaModel * do_clone() const override {
                return new RAWSOmegaModel(*this);
            }

            HOST_DEVICE_PREFIX void do_print() const override {
                printf("RAWSOmegaModel Control Model\n");
            }

            HOST_DEVICE_PREFIX void do_applyControl(Discretisation * const discretisation, const uint_t timeStep,
                                                    const real_t rotorMeanWind, const real_t rotorTorque,
                                                    const real_t rotorOmega, const Vector3<real_t> rotorWindVaneWind) override {
                previousRAWS_ = alphaFilter_ * previousRAWS_ + real_t(1-alphaFilter_)*rotorMeanWind;
                auto omegas = interpolator_(previousRAWS_);
                discretisation->setRelativeRotationalVelocity( omegas );
            }

            HOST_DEVICE_PREFIX void do_checkNeeds(bool & needsMeanVelocity, bool & needsTorque, bool & needsWindVane) override {
                needsMeanVelocity = true;
                needsTorque = false;
                needsWindVane = false;
            }

            HOST_DEVICE_PREFIX walberla::mpi::MPISize do_getItemSize() override {
                return 2 * sizeof(real_t);
            };

            HOST_DEVICE_PREFIX void do_fillBuffer(BufferType_T * buffer) override {

                memcpy(buffer, &alphaFilter_, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(buffer, &previousRAWS_, sizeof(real_t));
                buffer = buffer + sizeof(real_t);

            };

            HOST_DEVICE_PREFIX void do_unfillBuffer(BufferType_T * buffer) override {

                memcpy(&alphaFilter_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&previousRAWS_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);

            };

            data_interpolation::DataInterpolator<real_t, Vector3<real_t>, false> interpolator_{};
            real_t alphaFilter_{};
            real_t previousRAWS_{};

        };

    } // namespace control_model

} // namespace turbine_core

#endif //TURBINECORE_CONTROLMODEL_RAWSOmega_H
