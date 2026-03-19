
#pragma once

#ifndef TURBINECORE_DISCON_TORQUE_MODEL_H
#define TURBINECORE_DISCON_TORQUE_MODEL_H

#include "ControlModel.h"
#include "conversion/Conversion.h"

namespace turbine_core {

    namespace control_model {

        class DISCONTorqueModel final: public ControlModel {

            friend class ControlModel;

        public:

            static constexpr Type TYPE = Type::DISCON_TORQUE_CONTROL;

            HOST_DEVICE_PREFIX DISCONTorqueModel() {}

            HOST_DEVICE_PREFIX DISCONTorqueModel(
                const real_t & cutinSpeed,
                const real_t & region_2_startingspeed,
                const real_t & region_2_endingspeed,
                const real_t & rated_generator_speed,
                const real_t & cutinTorque,
                const real_t & ratedGeneratorTorque,
                const real_t & rotorRadius,
                const real_t & Ngear,
                const real_t & ratedPower,
                const real_t & CpOpt,
                const real_t & tsrOpt,
                const real_t & iRot,
                const real_t & iGen,
                const real_t & gearboxEfficiency,
                const real_t & generatorEfficiency
                ) {
                    cutinSpeed_ = cutinSpeed;
                    region_2_startingspeed_ = region_2_startingspeed;
                    region_2_endingspeed_ = region_2_endingspeed;
                    rated_generator_speed_ = rated_generator_speed;

                    cutinTorque_ = cutinTorque;
                    ratedGeneratorTorque_ = ratedGeneratorTorque;

                    rotorRadius_ = rotorRadius;
                    Ngear_ = Ngear;  // Gearbox ratio
                    ratedPower_ = ratedPower;
                    CpOpt_ = CpOpt;
                    tsrOpt_ = tsrOpt;
                    iRot_ = iRot;  // Rotor inertia
                    iGen_ = iGen; // Generator inertia
                    gearboxEfficiency_ = gearboxEfficiency;
                    generatorEfficiency_ = generatorEfficiency;
            }

            HOST_DEVICE_PREFIX DISCONTorqueModel( const DISCONTorqueModel & rawsModel)
                :
                  cutinSpeed_(rawsModel.cutinSpeed_),
                  region_2_startingspeed_(rawsModel.region_2_startingspeed_),
                  region_2_endingspeed_(rawsModel.region_2_endingspeed_),
                  rated_generator_speed_(rawsModel.rated_generator_speed_),
                  cutinTorque_(rawsModel.cutinTorque_),
                  ratedGeneratorTorque_(rawsModel.ratedGeneratorTorque_),
                  rotorRadius_(rawsModel.rotorRadius_),
                  Ngear_(rawsModel.Ngear_),
                  ratedPower_(rawsModel.ratedPower_),
                  CpOpt_(rawsModel.CpOpt_),
                  tsrOpt_(rawsModel.tsrOpt_),
                  iRot_(rawsModel.iRot_),
                  iGen_(rawsModel.iGen_),
                  gearboxEfficiency_(rawsModel.gearboxEfficiency_),
                  generatorEfficiency_(rawsModel.generatorEfficiency_)

            {}

            HOST_DEVICE_PREFIX DISCONTorqueModel & operator=( DISCONTorqueModel && model ) = delete;

            HOST_DEVICE_PREFIX ~DISCONTorqueModel() override = default;

            HOST_DEVICE_PREFIX auto & previousFilteredHighSpeedShaftSpeedSI() {
                return previousFilteredHighSpeedShaftSpeedSI_;
            }

            HOST_DEVICE_PREFIX auto & cutinSpeed() {
                return cutinSpeed_;
            }

            HOST_DEVICE_PREFIX auto & region_2_startingspeed() {
                return region_2_startingspeed_;
            }

            HOST_DEVICE_PREFIX auto & region_2_endingspeed() {
                return region_2_endingspeed_;
            }

            HOST_DEVICE_PREFIX auto & rated_generator_speed() {
                return rated_generator_speed_;
            }

            HOST_DEVICE_PREFIX auto & cutinTorque() {
                return cutinTorque_;
            }

            HOST_DEVICE_PREFIX auto & ratedGeneratorTorque() {
                return ratedGeneratorTorque_;
            }

            HOST_DEVICE_PREFIX auto & rotorRadius() {
                return rotorRadius_;
            }

            HOST_DEVICE_PREFIX auto & Ngear() {
                return Ngear_;
            }

            HOST_DEVICE_PREFIX auto & ratedPower() {
                return ratedPower_;
            }

            HOST_DEVICE_PREFIX auto & CpOpt() {
                return CpOpt_;
            }

            HOST_DEVICE_PREFIX auto & tsrOpt() {
                return tsrOpt_;
            }

            HOST_DEVICE_PREFIX auto & iRot() {
                return iRot_;
            }

            HOST_DEVICE_PREFIX auto & iGen() {
                return iGen_;
            }

            HOST_DEVICE_PREFIX auto & gearboxEfficiency() {
                return gearboxEfficiency_;
            }

            HOST_DEVICE_PREFIX auto & generatorEfficiency() {
                return generatorEfficiency_;
            }

        private:

            HOST_DEVICE_PREFIX real_t kOpt() {
                return (walberla::math::pi * real_t(1.2) * std::pow(rotorRadius_, real_t(5))
                    * CpOpt_) / (real_t(2) * (Ngear_*Ngear_*Ngear_) *
                    (tsrOpt_*tsrOpt_*tsrOpt_) * gearboxEfficiency_);
            }

            HOST_DEVICE_PREFIX real_t region_1_control() {
                return 0;
            }

            HOST_DEVICE_PREFIX real_t region_1half_control(real_t filteredHighSpeedShaftSpeed) {
                real_t d_generator_speed = filteredHighSpeedShaftSpeed - cutinSpeed_;
                real_t region_2_startingtorque = kOpt() * region_2_startingspeed_*region_2_startingspeed_;
                real_t torqueSlope = (region_2_startingtorque - cutinTorque_)/(region_2_startingspeed_ - cutinSpeed_);
                return cutinTorque_ + torqueSlope * d_generator_speed;
            }

            HOST_DEVICE_PREFIX real_t region_2_control(real_t filteredHighSpeedShaftSpeed) {
                return kOpt() * filteredHighSpeedShaftSpeed*filteredHighSpeedShaftSpeed;
            }

            HOST_DEVICE_PREFIX real_t region_2half_control(real_t filteredHighSpeedShaftSpeed) {
                real_t d_generator_speed = filteredHighSpeedShaftSpeed - region_2_endingspeed_;
                real_t region_2_endingtorque = kOpt() * region_2_endingspeed_*region_2_endingspeed_;
                real_t torque_slope = (ratedGeneratorTorque_ - region_2_endingtorque) /
                 (rated_generator_speed_ - region_2_endingspeed_);
                real_t generator_torque = region_2_endingtorque + torque_slope * d_generator_speed;
                return generator_torque;
            }

            HOST_DEVICE_PREFIX real_t region_3_control(real_t filteredHighSpeedShaftSpeed) {
                return ratedPower_ / filteredHighSpeedShaftSpeed;
            }

            HOST_DEVICE_PREFIX real_t calculateGeneratorTORQUE(real_t filteredHighSpeedShaftSpeed) {

                real_t generatorTorque;

                if ((cutinSpeed_ <= filteredHighSpeedShaftSpeed) && (filteredHighSpeedShaftSpeed <= region_2_startingspeed_)) {
                    generatorTorque = region_1half_control(filteredHighSpeedShaftSpeed);
                }
                else if ( (region_2_startingspeed_ <= filteredHighSpeedShaftSpeed) && (filteredHighSpeedShaftSpeed <= region_2_endingspeed_) ) {
                    generatorTorque = region_2_control(filteredHighSpeedShaftSpeed);
                }
                else if ( (region_2_endingspeed_ <= filteredHighSpeedShaftSpeed) &&  (filteredHighSpeedShaftSpeed <= rated_generator_speed_) ) {
                    generatorTorque = region_2half_control(filteredHighSpeedShaftSpeed);
                }
                else if (filteredHighSpeedShaftSpeed > rated_generator_speed_) {
                    generatorTorque = region_3_control(filteredHighSpeedShaftSpeed);
                }
                else {
                    generatorTorque = region_1_control();
                }

                return generatorTorque;
            }

            HOST_DEVICE_PREFIX DISCONTorqueModel * do_clone() const override {
                return new DISCONTorqueModel(*this);
            }

            HOST_DEVICE_PREFIX void do_print() const override {
                printf("DISCONTorqueModel Control Model\n");
            }

            HOST_DEVICE_PREFIX void do_applyControl(Discretisation * const discretisation, const uint_t timeStep,
                                                    const real_t rotorMeanWind, const real_t rotorTorque,
                                                    const real_t rotorOmega, const Vector3<real_t> rotorWindVaneWind) override {

                real_t filteredHighSpeedShaftSpeed = fabs(rotorOmega) * Ngear_;
                real_t generatorTorque = calculateGeneratorTORQUE(filteredHighSpeedShaftSpeed);

                // Per unit time step
                real_t domegaDt = (-rotorTorque * gearboxEfficiency_ - Ngear_ * generatorTorque) /
                    (iRot_ + Ngear_*Ngear_ * iGen_);
                filteredHighSpeedShaftSpeed += domegaDt * Ngear_;

                // Use low-speed shaft speed
                Vector3<real_t> omegas{0., 0., -real_t(filteredHighSpeedShaftSpeed / Ngear_ )};
                discretisation->setRelativeRotationalVelocity( omegas);
            }

            HOST_DEVICE_PREFIX void do_checkNeeds(bool & needsMeanVelocity, bool & needsTorque, bool & needsWindVane) override {
                needsMeanVelocity = false;
                needsTorque = true;
                needsWindVane = false;
            }

            HOST_DEVICE_PREFIX walberla::mpi::MPISize do_getItemSize() override {
                return 16 * sizeof(real_t);
            };

            HOST_DEVICE_PREFIX void do_fillBuffer(BufferType_T * buffer) override {
                //
                // fill buffer
                memcpy(buffer, &previousFilteredHighSpeedShaftSpeedSI_, sizeof(real_t));
                buffer = buffer + sizeof(previousFilteredHighSpeedShaftSpeedSI_);
                memcpy(buffer, &cutinSpeed_, sizeof(real_t));
                buffer = buffer + sizeof(cutinSpeed_);
                memcpy(buffer, &region_2_startingspeed_, sizeof(real_t));
                buffer = buffer + sizeof(region_2_startingspeed_);
                memcpy(buffer, &region_2_endingspeed_, sizeof(real_t));
                buffer = buffer + sizeof(region_2_endingspeed_);
                memcpy(buffer, &rated_generator_speed_, sizeof(real_t));
                buffer = buffer + sizeof(rated_generator_speed_);
                memcpy(buffer, &cutinTorque_, sizeof(real_t));
                buffer = buffer + sizeof(cutinTorque_);
                memcpy(buffer, &ratedGeneratorTorque_, sizeof(real_t));
                buffer = buffer + sizeof(ratedGeneratorTorque_);
                memcpy(buffer, &rotorRadius_, sizeof(real_t));
                buffer = buffer + sizeof(rotorRadius_);
                memcpy(buffer, &Ngear_, sizeof(real_t));
                buffer = buffer + sizeof(Ngear_);
                memcpy(buffer, &ratedPower_, sizeof(real_t));
                buffer = buffer + sizeof(ratedPower_);
                memcpy(buffer, &CpOpt_, sizeof(real_t));
                buffer = buffer + sizeof(CpOpt_);
                memcpy(buffer, &tsrOpt_, sizeof(real_t));
                buffer = buffer + sizeof(tsrOpt_);
                memcpy(buffer, &iRot_, sizeof(real_t));
                buffer = buffer + sizeof(iRot_);
                memcpy(buffer, &iGen_, sizeof(real_t));
                buffer = buffer + sizeof(iGen_);
                memcpy(buffer, &gearboxEfficiency_, sizeof(real_t));
                buffer = buffer + sizeof(gearboxEfficiency_);
                memcpy(buffer, &generatorEfficiency_, sizeof(real_t));
                buffer = buffer + sizeof(generatorEfficiency_);
            };

            HOST_DEVICE_PREFIX void do_unfillBuffer(BufferType_T * buffer) override {
                //
                // unfill buffer
                memcpy(&previousFilteredHighSpeedShaftSpeedSI_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&cutinSpeed_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&region_2_startingspeed_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&region_2_endingspeed_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&rated_generator_speed_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&cutinTorque_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&ratedGeneratorTorque_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&rotorRadius_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&Ngear_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&ratedPower_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&CpOpt_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&tsrOpt_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&iRot_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&iGen_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&gearboxEfficiency_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
                memcpy(&generatorEfficiency_, buffer, sizeof(real_t));
                buffer = buffer + sizeof(real_t);
            };

            real_t previousFilteredHighSpeedShaftSpeedSI_{real_t(.4*0.05)*real_t(96.)};
            real_t cutinSpeed_{};
            real_t region_2_startingspeed_{};
            real_t region_2_endingspeed_{};
            real_t rated_generator_speed_{};

            real_t cutinTorque_{};
            real_t ratedGeneratorTorque_{};

            real_t rotorRadius_{};
            real_t Ngear_{};  // Gearbox ratio
            real_t ratedPower_{};
            real_t CpOpt_{};
            real_t tsrOpt_{};
            real_t iRot_{};  // Rotor inertia
            real_t iGen_{}; // Generator inertia
            real_t gearboxEfficiency_{};
            real_t generatorEfficiency_{};
        };

    } // namespace control_model

} // namespace turbine_core

#endif //TURBINECORE_DISCON_TORQUE_MODEL_H
