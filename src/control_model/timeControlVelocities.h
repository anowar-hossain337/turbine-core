
#pragma once

#ifndef TURBINECORE_CONTROLMODEL_TimeControlVelocities_H
#define TURBINECORE_CONTROLMODEL_TimeControlVelocities_H

#include "ControlModel.h"
#include "conversion/Conversion.h"

namespace turbine_core {

    namespace control_model {

        class TimeControlVelocitiesModel final: public ControlModel {

            friend class ControlModel;

        public:

            static constexpr Type TYPE = Type::TIME_CONTROL_VELOCITIES;

            HOST_DEVICE_PREFIX TimeControlVelocitiesModel() = default;

            HOST_PREFIX explicit TimeControlVelocitiesModel(const std::string & controlTable) {

                std::ifstream file(controlTable);

                std::string line;

                std::vector<real_t> times;
                std::vector<Vector3<real_t>> translationVelocities, rotationVelocities;

                real_t C_u = Conversion::C_t() / Conversion::C_l();

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
                        real_t time, transX, transY, transZ, rotX, rotY, rotZ;
                        in >> time >> transX >> transY >> transZ >> rotX >> rotY >> rotZ;
                        WALBERLA_LOG_DETAIL_ON_ROOT("Phy LBM times " << time << ", " << time / Conversion::C_t())
                        times.push_back(time / Conversion::C_t());
                        translationVelocities.emplace_back(transX * C_u, transY * C_u, transZ * C_u);
                        real_t fix = 1./1.;
                        rotationVelocities.emplace_back(rotX * Conversion::C_t() * fix, rotY * Conversion::C_t()*fix, rotZ * Conversion::C_t()*fix);
                        WALBERLA_LOG_DETAIL_ON_ROOT("trans and rots " << transX * C_u << ", " << transY * C_u << ", " << transZ * C_u
                        << ", " << rotX * Conversion::C_t() << ", " << rotY * Conversion::C_t() << ", " << rotZ * Conversion::C_t())
                    }

                    file.close();
                }

                interpolatorTranslations_.setPoints(times.size(), times.data(), translationVelocities.data());
                interpolatorRotations_.setPoints(times.size(), times.data(), rotationVelocities.data());
            }

            HOST_DEVICE_PREFIX TimeControlVelocitiesModel( const uint_t nPoints,
                                               real_t * times, Vector3<real_t> * translationalVelocities,
                                               Vector3<real_t> * rotationalVelocities )
            {
                interpolatorTranslations_.setPoints(nPoints, times, translationalVelocities);
                interpolatorRotations_.setPoints(nPoints, times, rotationalVelocities);
            }

            HOST_DEVICE_PREFIX TimeControlVelocitiesModel( const TimeControlVelocitiesModel & model)
                    :  interpolatorTranslations_(model.interpolatorTranslations_),
                    interpolatorRotations_(model.interpolatorRotations_) {}

            HOST_DEVICE_PREFIX TimeControlVelocitiesModel & operator=( RAWSAnglesModel && model ) = delete;

            HOST_DEVICE_PREFIX ~TimeControlVelocitiesModel() override = default;

            HOST_DEVICE_PREFIX data_interpolation::DataInterpolator<real_t, Vector3<real_t>, false> * interpolatorTranslations() {
                return &interpolatorTranslations_;
            }

            HOST_DEVICE_PREFIX data_interpolation::DataInterpolator<real_t, Vector3<real_t>, false> * interpolatorRotations() {
                return &interpolatorRotations_;
            }

        private:

            HOST_DEVICE_PREFIX TimeControlVelocitiesModel * do_clone() const override {
                return new TimeControlVelocitiesModel(*this);
            }

            HOST_DEVICE_PREFIX void do_print() const override {
                printf("RAWSAnglesModel Control Model\n");
            }

            HOST_DEVICE_PREFIX void do_applyControl(Discretisation * const discretisation, const uint_t timeStep, 
                                                    const real_t rotorMeanWind, const real_t rotorTorque, 
                                                    const real_t rotorOmega, const Vector3<real_t> rotorWindVaneWind) override {

                auto translations = interpolatorTranslations_(timeStep);
//                printf("Set translations: %f, %f, %f: \n", translations[0], translations[1], translations[2]);
                auto rotations = interpolatorRotations_(timeStep);
                discretisation->setRelativeTranslationVelocity(translations);
                discretisation->setRelativeRotationalVelocity(rotations);
//                printf("Set rotations: %f, %f, %f: \n", rotations[0], rotations[1], rotations[2]);
            }

            HOST_DEVICE_PREFIX void do_checkNeeds(bool & needsMeanVelocity, bool & needsTorque, bool & needsWindVane) override {
                needsMeanVelocity = false;
                needsTorque = false;
                needsWindVane = false;
            };

            HOST_DEVICE_PREFIX walberla::mpi::MPISize do_getItemSize() override {return 0; };

            HOST_DEVICE_PREFIX void do_fillBuffer(BufferType_T * buffer) override {};

            HOST_DEVICE_PREFIX void do_unfillBuffer(BufferType_T * buffer) override {};

            data_interpolation::DataInterpolator<real_t, Vector3<real_t>, false> interpolatorTranslations_{};
            data_interpolation::DataInterpolator<real_t, Vector3<real_t>, false> interpolatorRotations_{};
        };

    } // namespace control_model

} // namespace turbine_core

#endif //TURBINECORE_CONTROLMODEL_TimeControlVelocities_H
