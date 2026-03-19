
#pragma once

#ifndef TURBINECORE_CONFIGTURBINECREATOR_H
#define TURBINECORE_CONFIGTURBINECREATOR_H

#include <cassert>
#include <utility>

#include <core/math/Vector3.h>
#include <core/math/Constants.h>
#include <core/logging/Logging.h>

#include "TurbineCreator.h"
#include "FileParser.h"

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/Enums.h"

#include "discretisation/all.h"
#include "force_model/all.h"
#include "control_model/all.h"

#include <core/config/Config.h>

namespace turbine_core {

    namespace creator {

        template< typename Topology_T, typename ForceField_T, typename Interpolator_T, typename ForceDistributor_T, typename TipLossModel_T >
        class ConfigTurbineCreator : public TurbineCreator<Topology_T> {

            using TopologyPtr = std::shared_ptr<Topology_T>;
            using ADM_T = force_model::ActuatorDiskModel<Interpolator_T, ForceDistributor_T>;

        public:

            HOST_PREFIX ConfigTurbineCreator(const std::string & outputDirectory_, const walberla::Config::BlockHandle & handle )
                    : blockHandle_(handle), outputDirectory_(outputDirectory_)
            {
                WALBERLA_LOG_INFO_ON_ROOT( "nTurbines = " << blockHandle_.getNumBlocks() )
            }

            HOST_PREFIX ~ConfigTurbineCreator() {}

        private:

            HOST_PREFIX void doCreateGeometry( std::vector<TopologyPtr> & turbines) const override {

                walberla::Config::Blocks turbineBlocks;
                blockHandle_.getBlocks(turbineBlocks);

                const auto nTurbines = turbineBlocks.size();
                turbines.reserve(nTurbines);
                for (auto & turbineBlock : turbineBlocks ) {

                    /// add basePoint
                    turbines.emplace_back(std::make_shared<Topology_T>());
                    auto& turbinePtr = turbines.back();
                    const std::string blockName = turbineBlock.getKey();

                    if(!outputDirectory_.empty()) {
                        WALBERLA_ROOT_SECTION() {
                            auto baseFolder = walberla::filesystem::path(outputDirectory_) / blockName;
                            // clear directory
                            if (walberla::filesystem::exists(baseFolder)) {
                                walberla::filesystem::remove_all(baseFolder);
                            }

                            if (!walberla::filesystem::exists(baseFolder)) {
                                walberla::filesystem::create_directories(baseFolder);
                            }
                        }
                    }

                    auto turbineConfig = blockHandle_.getBlock(blockName);
                    const Vector3<real_t> origin = turbineConfig.getParameter<Vector3<real_t>>("basePoint");
                    turbinePtr->setupEnvironment(blockName, origin);

                    const auto aeroModel = AeroModel::toType(turbineConfig.getParameter<std::string>( "aeroModel", "blades" ));

                    WALBERLA_LOG_PROGRESS_ON_ROOT("Adding geometry for turbine \" " << blockName << "\".")

                    addTowerGeometry(turbinePtr, turbineConfig);
                    addNacelleGeometry(turbinePtr, turbineConfig);
                    addHubGeometry(turbinePtr, turbineConfig);

                    if(aeroModel == AeroModel::ROTOR_DISK){
                        addRotorDiskGeometry(turbinePtr, turbineConfig);
                    } else {
                        addBladeGeometry(turbinePtr, turbineConfig);
                    }

                    turbinePtr->callback(Component::UPDATE_DISCRETISATION, ComponentType::ANY);
                }
            }

            HOST_PREFIX void doInitialiseForceAndControlModel( std::vector<TopologyPtr> & turbines, const std::vector<AABB> & turbineAABBs,
                                                               const std::shared_ptr<domain::TurbineDomain> & domain,
                                                               std::vector<bool> &isControlled,
                                                               std::vector <uint_t> &nPointsComponents,
                                                               std::vector <uint_t> &nComponents,
                                                               std::vector <uint_t> &nPointsComponentsForSpread) const override {

                walberla::Config::Blocks turbineBlocks;
                blockHandle_.getBlocks(turbineBlocks);

                isControlled.resize(turbines.size());
                nPointsComponents.resize(turbines.size());
                nComponents.resize(turbines.size());
                nPointsComponentsForSpread.resize(turbines.size());

                for (uint_t t = 0; t < turbines.size(); ++t) {

                    auto & localAABBs = domain->localAABBs();

                    bool containedLocally {false};

                    for (auto & localAABB : localAABBs) {
                        if(turbineAABBs[t].intersects(localAABB)) {
                            containedLocally = true;
                        }
                    }

                    bool turbineIsControlled{false};

                    if(!containedLocally) {
                        turbines[t].reset();
                    } else {

                        auto blockName = turbineBlocks[t].getKey();

                        /// create discretisation and force model
                        const auto turbineConfig = blockHandle_.getBlock(blockName);
                        const auto aeroModel = AeroModel::toType(turbineConfig.getParameter<std::string>( "aeroModel", "blades" ));

                        WALBERLA_LOG_PROGRESS_ON_ROOT("Adding force model for turbine \" " << blockName << "\" using " << aeroModel <<  ".")

                        /// TOWER
                        addTowerForceAndControlModel(turbines[t], turbineConfig, nPointsComponents[t], nComponents[t],
                                                     nPointsComponentsForSpread[t]);
                        addNacelleForceAndControlModel(turbines[t], turbineConfig, nPointsComponents[t],
                                                       nComponents[t]);
                        addHubForceAndControlModel(turbines[t], turbineConfig, nPointsComponents[t], nComponents[t],
                                                   nPointsComponentsForSpread[t]);

                        if(aeroModel == AeroModel::ROTOR_DISK) {
                            WALBERLA_CHECK(turbineConfig.getParameter<bool>("nacelleEffect") == false, "Nacelle effet cannot be combined with the rotor disk model.")
                            addRotorDiskForceAndControlModel(turbines[t], turbineConfig, nPointsComponents[t],
                                                             nComponents[t], nPointsComponentsForSpread[t]);
                        } else {
                            addBladeForceAndControlModel(turbines[t], turbineConfig, nPointsComponents[t],
                                                         nComponents[t], nPointsComponentsForSpread[t]);
                        }

                        if(turbines[t]->isControlled())
                            turbineIsControlled = true;

                    }

                    walberla::mpi::allReduceInplace(turbineIsControlled, walberla::mpi::LOGICAL_OR);

                    isControlled[t] = turbineIsControlled;

                }
            }

            HOST_PREFIX void addTowerGeometry(const TopologyPtr & turbine, const walberla::Config::BlockHandle & config ) const {

                const bool towerEffect = config.getParameter<bool>( "towerEffect", false );

                const uint_t nPoints = towerEffect ? config.getParameter<uint_t>("nPointsTower", uint_t(20)) : uint_t(0);
                const std::string towerFile = config.getParameter<std::string>( "towerDescription" );

                WALBERLA_LOG_INFO_ON_ROOT("TowerEffect = " << towerEffect << ", nPoints = " << nPoints )

                auto towerPositions = std::make_shared<Vector3DataInterpolator<false>>();
                auto towerAngles = std::make_shared<Vector3DataInterpolator<false>>();
                auto towerWidths = std::make_shared<ScalarDataInterpolator<false>>();
                std::vector<aerodynamics::AirfoilPolar> airfoilPolars;

                parseDescriptionFile(towerFile, nPoints, towerPositions, towerAngles, towerWidths, airfoilPolars);

                // set relative position/ velocity/ etc of discretisation
                const auto towerOrientation = Quaternion<real_t>( Vector3<real_t>{0,0,1}, walberla::math::pi );
                const Vector3<real_t> towerRotation{};
                const Vector3<real_t> towerTranslation{};
//                Vector3<real_t> towerTranslation(0,-0.05,0);

                auto discretisation = std::make_shared<discretisation::Line>(Vector3<real_t>(), towerTranslation, towerRotation, towerOrientation);
                discretisation->setPoints(nPoints, towerPositions.get(), towerAngles.get());

                turbine->addDiscretisation(config.getKey(), "Tower", ComponentType::TOWER, discretisation);

            } // addTowerGeometry

            HOST_PREFIX void
            addTowerForceAndControlModel(const TopologyPtr &turbine, const walberla::Config::BlockHandle &config,
                                         uint_t &nPointsComponents,
                                         uint_t &nComponents, uint_t &nPointsComponentsForSpread) const {

                const bool towerEffect = config.getParameter<bool>( "towerEffect", false );
//                WALBERLA_ASSERT(!towerEffect, "towerEffect is not supported")

                const uint_t nPoints = towerEffect ? config.getParameter<uint_t>("nPointsTower", uint_t(20)) : uint_t(0);
                nPointsComponents += nPoints;
                nComponents += 1;
                nPointsComponentsForSpread += nPoints;

                const std::string towerFile = config.getParameter<std::string>( "towerDescription" );

                WALBERLA_LOG_PROGRESS_ON_ROOT("TowerEffect = " << towerEffect << ", nPoints = " << nPoints )

                auto towerPositions = std::make_shared<Vector3DataInterpolator<false>>();
                auto towerAngles = std::make_shared<Vector3DataInterpolator<false>>();
                auto towerWidths = std::make_shared<ScalarDataInterpolator<false>>();
                std::vector<aerodynamics::AirfoilPolar> airfoilPolars;

                parseDescriptionFile(towerFile, nPoints, towerPositions, towerAngles, towerWidths, airfoilPolars);

                auto nInterpolationPoints = towerPositions->length();

                const auto towerControl = control_model::ControlModel::toType(
                    config.getParameter<std::string>( "towerController", "None" ));

                if(!towerEffect) {

                    std::shared_ptr<ForceModel> forceModel = nullptr;

                    if(towerControl == control_model::ControlModel::TIME_CONTROL_VELOCITIES) {
                        auto towerControlBlock = config.getOneBlock(control_model::ControlModel::toString(towerControl));
                        const std::string controlTable = towerControlBlock.template getParameter<std::string>("table");
                        auto controlModel = std::make_shared<control_model::TimeControlVelocitiesModel>(controlTable);

                        turbine->initialiseForceAndControlModel("Tower", forceModel, controlModel);
                    } else {
                        std::shared_ptr<control_model::ControlModel> controlModel = nullptr;
                        turbine->initialiseForceAndControlModel("Tower", forceModel, controlModel);
                    }

                } else {

                    force_model::NoTipLoss tipLossModel{};
                    auto forceModel = std::make_shared<force_model::ActuatorLineModel<Interpolator_T, ForceDistributor_T, force_model::NoTipLoss>>(nPoints, nInterpolationPoints, towerPositions->x(), airfoilPolars.data(), *towerWidths, tipLossModel);
                    
                    const auto towerControl = control_model::ControlModel::toType(
                        config.getParameter<std::string>( "towerController", "None" ));
                    if(towerControl == control_model::ControlModel::TIME_CONTROL_VELOCITIES) {
                        auto towerControlBlock = config.getOneBlock(control_model::ControlModel::toString(towerControl));
                        const std::string controlTable = towerControlBlock.template getParameter<std::string>("table");
                        auto controlModel = std::make_shared<control_model::TimeControlVelocitiesModel>(controlTable);
                        
                        turbine->initialiseForceAndControlModel("Tower", forceModel, controlModel);
                    } else {
                        std::shared_ptr<control_model::ControlModel> controlModel = nullptr;
                        
                        turbine->initialiseForceAndControlModel("Tower", forceModel, controlModel);
                    }

                } // else (!towerEffect)

            } // addTowerForceAndControlModel

            HOST_PREFIX void addNacelleGeometry(const TopologyPtr & turbine, const walberla::Config::BlockHandle & config ) const {

                const real_t nacelleYaw = config.getParameter<real_t>("nacelleYaw");

                const real_t hubRadius = real_t(config.getParameter<real_t>( "hubRadius_SI" )) / Conversion::C_l();

                Vector3<real_t> relPosNacelle (0,0,hubRadius);
                Quaternion<real_t> relOrientNacelle = Quaternion<real_t>(Vector3<real_t>(1,0,0), nacelleYaw * walberla::math::pi / real_t(180));
                relOrientNacelle = Quaternion<real_t>(Vector3<real_t>(0,1,0), walberla::math::half_pi) * relOrientNacelle;

                auto discretisation = std::make_shared<discretisation::Disk>(hubRadius, relPosNacelle, Vector3<real_t>(), relOrientNacelle);

                turbine->addDiscretisation("Tower", "Nacelle", ComponentType::NACELLE, discretisation);

            }

            HOST_PREFIX void
            addNacelleForceAndControlModel(const TopologyPtr &turbine, const walberla::Config::BlockHandle &config,
                                           uint_t &nPointsComponents,
                                           uint_t &nComponents) const {

                nPointsComponents += 1;
                nComponents += 1;
                std::shared_ptr<ForceModel> forceModel = nullptr;

                const auto nacelleControl = control_model::ControlModel::toType(
                    config.getParameter<std::string>("nacelleController", "None"));

                    if (nacelleControl == control_model::ControlModel::SIMPLE_YAW_CONTROLLER) {
                        auto nacelleControlBlock = config.getOneBlock(control_model::ControlModel::toString(nacelleControl));
                        real_t yawingSpeed = nacelleControlBlock.template getParameter<real_t>("yawingSpeed");
                        auto controlModel = std::make_shared<control_model::SimpleYawController>(yawingSpeed);
                        turbine->initialiseForceAndControlModel("Nacelle", forceModel, controlModel);
                    } else {
                        std::shared_ptr<control_model::ControlModel> controlModel = nullptr;
                        turbine->initialiseForceAndControlModel("Nacelle", forceModel, controlModel);
                    }
            }

            HOST_PREFIX void addRotorDiskGeometry(const TopologyPtr & turbine, const walberla::Config::BlockHandle & config) const {

                const real_t rotorRadius = config.getParameter<real_t>("diameter_LU") / real_t(2);

                Vector3<real_t> relPosRotorDisk (0,0,0);
                Quaternion<real_t> relOrientRotorDisk = Quaternion<real_t>(Vector3<real_t>(0,1,0), 0. * walberla::math::pi / real_t(180));

                auto discretisation = std::make_shared<discretisation::Disk>(rotorRadius, relPosRotorDisk, Vector3<real_t>(), relOrientRotorDisk);
                std::shared_ptr<control_model::ControlModel> controlModel = nullptr;

                turbine->addDiscretisation("Nacelle", "RotorDisk", ComponentType::ROTOR_DISK, discretisation);

            }

            HOST_PREFIX void
            addRotorDiskForceAndControlModel(const TopologyPtr &turbine, const walberla::Config::BlockHandle &config,
                                             uint_t &nPointsComponents,
                                             uint_t &nComponents, uint_t &nPointsComponentsForSpread) const {

                const real_t rotorDiameter = config.getParameter<real_t>("diameter_LU");
                nComponents += 1;
                nPointsComponents += 1;
                nPointsComponentsForSpread += std::pow(rotorDiameter, 3);


                real_t rotorDrag = config.getParameter<real_t>("rotorDrag");
                auto forceModel = std::make_shared<ADM_T>(rotorDiameter / real_t(2), rotorDrag);
                std::shared_ptr<control_model::ControlModel> controlModel = nullptr;

                turbine->initialiseForceAndControlModel("RotorDisk", forceModel, controlModel);
            }

            HOST_PREFIX void addHubGeometry( const TopologyPtr & turbine, const walberla::Config::BlockHandle & config ) const {

                const real_t hubRadius = real_t(config.getParameter<real_t>( "hubRadius_SI" )) / Conversion::C_l();
                const real_t hubDeport = real_t(config.getParameter<real_t>( "hubDeport_SI" )) / Conversion::C_l();

                Vector3<real_t> relPosHub (0,0,hubDeport);
                Vector3<real_t> hubRotation = config.getParameter<Vector3<real_t>>("hubRotationalVelocity");
                hubRotation = hubRotation * Conversion::C_t();

                hubRotation *= real_t(-1.0);

                const real_t hubTilt = config.getParameter<real_t>( "hubTilt" );
                Quaternion<real_t> relOrientHub = Quaternion<real_t>(Vector3<real_t>(0,1,0), -hubTilt * walberla::math::pi / real_t(180));

                auto discretisation = std::make_shared<discretisation::Disk>(hubRadius, relPosHub, hubRotation, relOrientHub);

                turbine->addDiscretisation("Nacelle", "Hub", ComponentType::HUB, discretisation);

            }

            HOST_PREFIX void
            addHubForceAndControlModel(const TopologyPtr &turbine, const walberla::Config::BlockHandle &config,
                                       uint_t &nPointsComponents,
                                       uint_t &nComponents, uint_t &nPointsComponentsForSpread) const {

                nPointsComponents += 1;
                nComponents += 1;

                const bool useNacelleEffect = config.getParameter<bool>("nacelleEffect", false);

                const real_t hubRadius = real_t(config.getParameter<real_t>("hubRadius_SI")) / Conversion::C_l();

                const auto hubControl = control_model::ControlModel::toType(
                    config.getParameter<std::string>("hubController", "None"));
                // HUB
                if (useNacelleEffect) {
                    real_t hubDrag{4.0};
                    auto forceModel = std::make_shared<ADM_T>(hubRadius, hubDrag);
                    nPointsComponentsForSpread += 8 * (std::pow(std::ceil(hubRadius), 3));


                    if (hubControl == control_model::ControlModel::RAWS_OMEGA) {
                        auto hubControlBlock = config.getOneBlock(control_model::ControlModel::toString(hubControl));
                        const std::string controlTable = hubControlBlock.template getParameter<std::string>("table");
                        real_t alphaFilter = hubControlBlock.template getParameter<real_t>("alphaFilter");
                        auto controlModel = std::make_shared<control_model::RAWSOmegaModel>(controlTable, alphaFilter);

                        turbine->initialiseForceAndControlModel("Hub", forceModel, controlModel);
                    } else if (hubControl == control_model::ControlModel::DISCON_TORQUE_CONTROL){
                        auto hubControlBlock = config.getOneBlock(control_model::ControlModel::toString(hubControl));

                        real_t C_torque = real_t(1) / (Conversion::C_m() * Conversion::C_l() * Conversion::C_l() / Conversion::C_t() / Conversion::C_t());
                        real_t C_inertia = real_t(1) / Conversion::C_m() / Conversion::C_l() / Conversion::C_l();
                        real_t cutinSpeed = hubControlBlock.template getParameter<real_t>("cutinSpeed") * Conversion::C_t();
                        real_t region_2_startingspeed = hubControlBlock.template getParameter<real_t>("region_2_startingspeed") * Conversion::C_t();
                        real_t region_2_endingspeed = hubControlBlock.template getParameter<real_t>("region_2_endingspeed") * Conversion::C_t();
                        real_t rated_generator_speed = hubControlBlock.template getParameter<real_t>("rated_generator_speed") * Conversion::C_t();
                        real_t cutin_torque = hubControlBlock.template getParameter<real_t>("cutin_torque") * C_torque;
                        real_t rated_generator_torque = hubControlBlock.template getParameter<real_t>("rated_generator_torque") * C_torque;
                        real_t rotor_radius = hubControlBlock.template getParameter<real_t>("rotor_radius") / Conversion::C_l();
                        real_t Ngear = hubControlBlock.template getParameter<real_t>("Ngear");
                        real_t rated_power = hubControlBlock.template getParameter<real_t>("rated_power") * C_torque * Conversion::C_t();
                        real_t Cp_opt = hubControlBlock.template getParameter<real_t>("Cp_opt");
                        real_t tsrOpt = hubControlBlock.template getParameter<real_t>("tsrOpt");
                        real_t Irot = hubControlBlock.template getParameter<real_t>("Irot") * C_inertia;
                        real_t Igen = hubControlBlock.template getParameter<real_t>("Igen") * C_inertia;
                        real_t gearboxEfficiency = hubControlBlock.template getParameter<real_t>("gearboxEfficiency");
                        real_t generatorEfficiency = hubControlBlock.template getParameter<real_t>("generatorEfficiency");

                        auto controlModel = std::make_shared<control_model::DISCONTorqueModel>(
                            cutinSpeed, region_2_startingspeed, region_2_endingspeed, rated_generator_speed,
                            cutin_torque, rated_generator_torque, rotor_radius, Ngear, rated_power, Cp_opt,
                            tsrOpt, Irot, Igen, gearboxEfficiency, generatorEfficiency
                            );

                        turbine->initialiseForceAndControlModel("Hub", forceModel, controlModel);
                    } else {
                        std::shared_ptr<control_model::ControlModel> controlModel = nullptr;

                        turbine->initialiseForceAndControlModel("Hub", forceModel, controlModel);
                    }
                } else {
                    std::shared_ptr<force_model::ForceModel> forceModel = nullptr;
                    if (hubControl == control_model::ControlModel::RAWS_OMEGA) {
                        auto hubControlBlock = config.getOneBlock(control_model::ControlModel::toString(hubControl));
                        const std::string controlTable = hubControlBlock.template getParameter<std::string>("table");
                        real_t alphaFilter = hubControlBlock.template getParameter<real_t>("alphaFilter");
                        auto controlModel = std::make_shared<control_model::RAWSOmegaModel>(controlTable, alphaFilter);

                        turbine->initialiseForceAndControlModel("Hub", forceModel, controlModel);
                    } else if (hubControl == control_model::ControlModel::DISCON_TORQUE_CONTROL){
                        auto hubControlBlock = config.getOneBlock(control_model::ControlModel::toString(hubControl));

                        real_t C_torque = real_t(1) / (Conversion::C_m() * Conversion::C_l() * Conversion::C_l() / Conversion::C_t() / Conversion::C_t());
                        real_t C_inertia = real_t(1) / Conversion::C_m() / Conversion::C_l() / Conversion::C_l();
                        real_t cutinSpeed = hubControlBlock.template getParameter<real_t>("cutinSpeed") * Conversion::C_t();
                        real_t region_2_startingspeed = hubControlBlock.template getParameter<real_t>("region_2_startingspeed") * Conversion::C_t();
                        real_t region_2_endingspeed = hubControlBlock.template getParameter<real_t>("region_2_endingspeed") * Conversion::C_t();
                        real_t rated_generator_speed = hubControlBlock.template getParameter<real_t>("rated_generator_speed") * Conversion::C_t();
                        real_t cutin_torque = hubControlBlock.template getParameter<real_t>("cutin_torque") * C_torque;
                        real_t rated_generator_torque = hubControlBlock.template getParameter<real_t>("rated_generator_torque") * C_torque;
                        real_t rotor_radius = hubControlBlock.template getParameter<real_t>("rotor_radius") / Conversion::C_l();
                        real_t Ngear = hubControlBlock.template getParameter<real_t>("Ngear");
                        real_t rated_power = hubControlBlock.template getParameter<real_t>("rated_power") * C_torque * Conversion::C_t();
                        real_t Cp_opt = hubControlBlock.template getParameter<real_t>("Cp_opt");
                        real_t tsrOpt = hubControlBlock.template getParameter<real_t>("tsrOpt");
                        real_t Irot = hubControlBlock.template getParameter<real_t>("Irot") * C_inertia;
                        real_t Igen = hubControlBlock.template getParameter<real_t>("Igen") * C_inertia;
                        real_t gearboxEfficiency = hubControlBlock.template getParameter<real_t>("gearboxEfficiency");
                        real_t generatorEfficiency = hubControlBlock.template getParameter<real_t>("generatorEfficiency");

                        auto controlModel = std::make_shared<control_model::DISCONTorqueModel>(
                            cutinSpeed, region_2_startingspeed, region_2_endingspeed, rated_generator_speed,
                            cutin_torque, rated_generator_torque, rotor_radius, Ngear, rated_power, Cp_opt,
                            tsrOpt, Irot, Igen, gearboxEfficiency, generatorEfficiency
                            );

                        turbine->initialiseForceAndControlModel("Hub", forceModel, controlModel);
                    } else {
                        std::shared_ptr<control_model::ControlModel> controlModel = nullptr;

                        turbine->initialiseForceAndControlModel("Hub", forceModel, controlModel);
                    }
                }
            }

            HOST_PREFIX void
            addBladeGeometry(const TopologyPtr &turbine, const walberla::Config::BlockHandle &config) const {

                const real_t toRad = walberla::math::pi / real_t(180.0);
                const real_t hubRadius = real_t(config.getParameter<real_t>( "hubRadius_SI" )) / Conversion::C_l();

                const uint_t nBlades = config.getParameter<uint_t>("nBlades");
                const std::string bladeFile = config.getParameter<std::string>("bladeDescription");

                const real_t bladePrecone = config.getParameter<real_t>("bladePrecone") * toRad;
                const real_t bladePitch = config.getParameter<real_t>("bladePitch") * toRad;

                Vector3<real_t> relPosBlade(-hubRadius,0,0);

                const uint_t nPointsBlades = config.getParameter<uint_t>("nPointsBlades", uint_t(40));


                Quaternion<real_t> pitchMatrix(Vector3<real_t>(0,0,1), -bladePitch);
                Quaternion<real_t> preconeMatrix(Vector3<real_t>(0,1,0), bladePrecone);

                auto halfPiY = Quaternion<real_t>(Vector3<real_t>(0,1,0), -walberla::math::half_pi);
                auto piZ = Quaternion<real_t>(Vector3<real_t>(0,0,1), walberla::math::pi);

                real_t azimuthOffset{ real_t(2.0) * walberla::math::pi / real_t(nBlades) };

                auto bladePositions = std::make_shared<Vector3DataInterpolator<false>>();
                auto bladeAngles = std::make_shared<Vector3DataInterpolator<false>>();
                auto bladeWidths = std::make_shared<ScalarDataInterpolator<false>>();
                std::vector<aerodynamics::AirfoilPolar> airfoilPolars;

                parseDescriptionFile(bladeFile, nPointsBlades, bladePositions, bladeAngles, bladeWidths, airfoilPolars);

                // set relative position/ velocity/ etc of discretisation

                for(uint_t b = 0; b < nBlades; ++b) {

                    real_t azimuth = azimuthOffset * real_t(b);

                    Quaternion<real_t> relOrientBlade = halfPiY * Quaternion<real_t>( Vector3<real_t>(1,0,0), -azimuth ) * preconeMatrix * pitchMatrix * piZ;

                    auto bladePosition = Quaternion<real_t>(Vector3<real_t>(0,0,1), -azimuth ).rotate(relPosBlade);

                    auto discretisation = std::make_shared<discretisation::Line>(bladePosition, Vector3<real_t>(), Vector3<real_t>(), relOrientBlade);
                    discretisation->setPoints(nPointsBlades, bladePositions.get(), bladeAngles.get());

                    if(!outputDirectory_.empty()) {

                        auto bladeFolder = walberla::filesystem::path(outputDirectory_).append(config.getKey()).append("Blade_" + std::to_string(b));

                        WALBERLA_ROOT_SECTION() {
                            // clear directory
                            if (walberla::filesystem::exists(bladeFolder)) {
                                walberla::filesystem::remove_all(bladeFolder);
                            }

                            if (!walberla::filesystem::exists(bladeFolder)) {
                                walberla::filesystem::create_directories(bladeFolder);
                            }

                            auto line = dynamic_cast<discretisation::Line*>(discretisation.get());
                            if(!line)
                                WALBERLA_ABORT("Blade must be of Line type.")

                            for (uint_t p = 0; p < nPointsBlades; ++p) {

                                const std::string pString{std::to_string(p)};

                                const std::string filename{ "Point_" + std::string(3 - pString.length(), '0') + pString + ".txt"};

                                std::ofstream os(bladeFolder / filename, std::ios::app);

                                os << "# FORCE OUTPUT for Point " << p << " "
                                   << "( spanwise location: "
                                   << ((line->points()[p].relativePosition).length() + line->relativePosition_.length()) * Conversion::C_l() << " )\n\n";

                                os << "# Timestep\tAzimuth angle\tNormal force\tTangential force\tAngle-of-attack\tLift\tDrag\tRelative velocity\n\n";

                                os.close();
                            }
                        }

                    }

                    turbine->addDiscretisation("Hub", "Blade_" + std::to_string(b), ComponentType::BLADE, discretisation);
                }

            } // add BladeGeometry

            HOST_PREFIX void
            addBladeForceAndControlModel(const TopologyPtr &turbine, const walberla::Config::BlockHandle &config,
                                         uint_t &nPointsComponents,
                                         uint_t &nComponents, uint_t &nPointsComponentsForSpread) const {

                const real_t toRad = walberla::math::pi / real_t(180.0);
                const real_t hubRadius = real_t(config.getParameter<real_t>("hubRadius_SI")) / Conversion::C_l();

                uint_t nBlades = config.getParameter<uint_t>("nBlades");
                nComponents += nBlades;


                const std::string bladeFile = config.getParameter<std::string>("bladeDescription");

                const real_t bladePrecone = config.getParameter<real_t>("bladePrecone") * toRad;
                const real_t bladePitch = config.getParameter<real_t>("bladePitch") * toRad;

                Vector3 <real_t> relPosBlade(-hubRadius, 0, 0);

                uint_t nPointsBlades = config.getParameter<uint_t>("nPointsBlades", uint_t(40));
                nPointsComponents += nPointsBlades;
                nPointsComponentsForSpread += nPointsBlades;


                Quaternion <real_t> pitchMatrix(Vector3<real_t>(0, 0, 1), -bladePitch);
                Quaternion <real_t> preconeMatrix(Vector3<real_t>(0, 1, 0), bladePrecone);

                auto halfPiY = Quaternion<real_t>(Vector3<real_t>(0, 1, 0), -walberla::math::half_pi);
                auto piZ = Quaternion<real_t>(Vector3<real_t>(0, 0, 1), walberla::math::pi);

                real_t azimuthOffset{real_t(2.0) * walberla::math::pi / real_t(nBlades)};

                auto bladePositions = std::make_shared<Vector3DataInterpolator < false>>();
                auto bladeAngles = std::make_shared<Vector3DataInterpolator < false>>();
                auto bladeWidths = std::make_shared<ScalarDataInterpolator < false>>();
                std::vector<aerodynamics::AirfoilPolar> airfoilPolars;

                parseDescriptionFile(bladeFile, nPointsBlades, bladePositions, bladeAngles, bladeWidths, airfoilPolars);

                auto nInterpolationPoints = bladePositions->length();

                // set relative position/ velocity/ etc of discretisation

                for (uint_t b = 0; b < nBlades; ++b) {

                    real_t azimuth = azimuthOffset * real_t(b);

                    Quaternion <real_t> relOrientBlade =
                            halfPiY * Quaternion<real_t>(Vector3<real_t>(1, 0, 0), -azimuth) * preconeMatrix *
                            pitchMatrix * piZ;

                    auto bladePosition = Quaternion<real_t>(Vector3<real_t>(0, 0, 1), -azimuth).rotate(
                            relPosBlade);

                    auto discretisation = std::make_shared<discretisation::Line>(bladePosition,
                                                                                 Vector3<real_t>(),
                                                                                 Vector3<real_t>(),
                                                                                 relOrientBlade);
                    discretisation->setPoints(nPointsBlades, bladePositions.get(), bladeAngles.get());


                    const bool useTipLoss = config.getParameter<bool>("tipLoss");

                    const auto bladeControl = control_model::ControlModel::toType(
                        config.getParameter<std::string>("bladeController", "None"));

                    if (bladeControl == control_model::ControlModel::RAWS_ANGLES) {
                        auto bladeControlBlock = config.getOneBlock(control_model::ControlModel::toString(bladeControl));
                        const std::string controlTable = bladeControlBlock.getParameter<std::string>("table");
                        const Vector3 <real_t> initialAngles = bladeControlBlock.getParameter < Vector3 < real_t >> ("initialAngles");
                        real_t alphaFilter = bladeControlBlock.template getParameter<real_t>("alphaFilter");
                        auto controlModel = std::make_shared<control_model::RAWSAnglesModel>(controlTable,
                                                                                             initialAngles,
                                                                                             alphaFilter);
                        if (useTipLoss) {
                            force_model::PrandtlTipLoss tipLossModel{};
                            auto forceModel = std::make_shared<force_model::ActuatorLineModel <Interpolator_T, ForceDistributor_T, force_model::PrandtlTipLoss>>
                            (nPointsBlades, nInterpolationPoints,
                                    bladePositions->x(), airfoilPolars.data(),
                                    *bladeWidths, tipLossModel);
                            turbine->initialiseForceAndControlModel("Blade_" + std::to_string(b), forceModel, controlModel);
                        } else {
                            force_model::NoTipLoss tipLossModel{};
                            auto forceModel = std::make_shared<force_model::ActuatorLineModel <Interpolator_T, ForceDistributor_T, force_model::NoTipLoss>>
                            (nPointsBlades, nInterpolationPoints,
                                    bladePositions->x(), airfoilPolars.data(),
                                    *bladeWidths, tipLossModel);
                            turbine->initialiseForceAndControlModel("Blade_" + std::to_string(b), forceModel, controlModel);
                        }
                    } else if (bladeControl == control_model::ControlModel::DISCON_PITCH_CONTROL) {
                        auto bladeControlBlock = config.getOneBlock(control_model::ControlModel::toString(bladeControl));
                        real_t pitchK = bladeControlBlock.template getParameter<real_t>("pitchK");
                        real_t pitchControlKP = bladeControlBlock.template getParameter<real_t>("pitchControlKP") / Conversion::C_t();
                        real_t pitchControlKI = bladeControlBlock.template getParameter<real_t>("pitchControlKI");
                        real_t pitchMin = bladeControlBlock.template getParameter<real_t>("pitchMin");
                        real_t pitchMax = bladeControlBlock.template getParameter<real_t>("pitchMax");
                        real_t maxPitchRate = bladeControlBlock.template getParameter<real_t>("maxPitchRate") * Conversion::C_t();
                        real_t ratedGeneratorSpeed = bladeControlBlock.template getParameter<real_t>("ratedGeneratorSpeed") * Conversion::C_t();
                        real_t nGear = bladeControlBlock.template getParameter<real_t>("nGear");

                        auto controlModel = std::make_shared<control_model::DISCONPitchModel>(
                                pitchK, pitchControlKP, pitchControlKI, pitchMin, pitchMax,
                                ratedGeneratorSpeed, nGear, maxPitchRate);

                        if (useTipLoss) {
                            force_model::PrandtlTipLoss tipLossModel{};
                            auto forceModel = std::make_shared<force_model::ActuatorLineModel <Interpolator_T, ForceDistributor_T, force_model::PrandtlTipLoss>>
                            (nPointsBlades, nInterpolationPoints,
                                    bladePositions->x(), airfoilPolars.data(),
                                    *bladeWidths, tipLossModel);
                            turbine->initialiseForceAndControlModel("Blade_" + std::to_string(b), forceModel, controlModel);
                        } else {
                            force_model::NoTipLoss tipLossModel{};
                            auto forceModel = std::make_shared<force_model::ActuatorLineModel <Interpolator_T, ForceDistributor_T, force_model::NoTipLoss>>
                            (nPointsBlades, nInterpolationPoints,
                                    bladePositions->x(), airfoilPolars.data(),
                                    *bladeWidths, tipLossModel);
                            turbine->initialiseForceAndControlModel("Blade_" + std::to_string(b), forceModel, controlModel);
                        }
                    } else {
                        std::shared_ptr<control_model::ControlModel> controlModel = nullptr;
                        if (useTipLoss) {
                            force_model::PrandtlTipLoss tipLossModel{};
                            auto forceModel = std::make_shared<force_model::ActuatorLineModel <Interpolator_T, ForceDistributor_T, force_model::PrandtlTipLoss>>
                                                (nPointsBlades, nInterpolationPoints,
                                                        bladePositions->x(), airfoilPolars.data(),
                                                        *bladeWidths, tipLossModel);
                            turbine->initialiseForceAndControlModel("Blade_" + std::to_string(b), forceModel, controlModel);
                        } else {
                            force_model::NoTipLoss tipLossModel{};
                            auto forceModel = std::make_shared<force_model::ActuatorLineModel <
                                    Interpolator_T, ForceDistributor_T, force_model::NoTipLoss>>
                            (nPointsBlades, nInterpolationPoints,
                                    bladePositions->x(), airfoilPolars.data(),
                                    *bladeWidths, tipLossModel);
                            turbine->initialiseForceAndControlModel("Blade_" + std::to_string(b), forceModel, controlModel);
                        }
                    }

                } // add Blades
            }

            const walberla::Config::BlockHandle blockHandle_;
            const std::string outputDirectory_;


        };

    }
}

#endif //TURBINECORE_CONFIGTURBINECREATOR_H
