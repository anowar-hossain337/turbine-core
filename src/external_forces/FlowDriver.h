//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file FlowDriver.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef TURBINECORE_FLOWDRIVER_H
#define TURBINECORE_FLOWDRIVER_H

#include <functional>

#include <core/Set.h>
#include <core/selectable/IsSetSelected.h>

#include <core/uid/SUID.h>

#include <domain_decomposition/IBlock.h>
#include <blockforest/BlockForest.h>
#include <blockforest/StructuredBlockForest.h>

#include "wind_turbine_core/WalberlaDataTypes.h"

#include "conversion/Conversion.h"

#include "data_interpolator/DataInterpolator.h"

#include "external_forces/DriverSetup.h"
#include "external_forces/ZeroForce.h"
#include "external_forces/PressureGradient.h"
#include "external_forces/DynamicPressureGradient.h"
#include "external_forces/CoriolisForce.h"

namespace turbine_core {

    namespace external_forces {

        template< typename T, typename Value_T >
        struct FlowDriverBase {

            FlowDriverBase( const std::shared_ptr<walberla::blockforest::StructuredBlockForest> & forest,
                            std::shared_ptr<walberla::Config> & globalConfig,
                            const domain::DomainSetup & domainSetup,
                            const BlockDataID forceFieldID, const BlockDataID velocityFieldID,
                            walberla::timeloop::ITimeloop * timeloop,
                            const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                            const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
                            : forest_(forest), forceFieldId_(forceFieldID), velocityFieldId_(velocityFieldID),
                              timeloop_(timeloop),
                              requiredBlockSelectors_(requiredBlockSelectors),
                              incompatibleBlockSelectors_(incompatibleBlockSelectors)
            {
                WALBERLA_LOG_WARNING_ON_ROOT("The manually implemented flow driver is deprecated and will be removed soon."
                                             "Please use the generated version implemented in python/generate_flow_driver.py")
                WALBERLA_CHECK(forest->storesUniformBlockGrid(), "The driving forces currently only work for uniform grids.")

                const auto defaultValue = std::numeric_limits<Value_T>::quiet_NaN();

                const auto parameters = globalConfig->getOneBlock("Parameters");
                const auto initialisation = globalConfig->getOneBlock("Initialisation");

                if(!domainSetup.periodicity_[0]) {
                    driverType_ = DriverSetup::None;
                    if(parameters.isDefined("driverType")) {
                        const auto type = DriverSetup::toType(parameters.getParameter<std::string>("driverType"));
                        if(type != DriverSetup::None)
                            WALBERLA_LOG_WARNING_ON_ROOT("Do you really want to set a flow driver for domains that are not periodic in the flow direction ?")
                        driverType_ = type;
                    }

                } else {
                    driverType_ = DriverSetup::toType(parameters.getParameter<std::string>("driverType"));
                }

                // Pressure-driven
                const Value_T kappa = parameters.getParameter<Value_T>("kappa", real_t(0.42));
                const Vector3<Value_T> initialVelocity = initialisation.getParameter<Vector3<Value_T>>("initialVelocity", Vector3<Value_T>(defaultValue));
                const real_t roughnessLengthRatio = parameters.getParameter<Value_T>("roughnessLengthRatio", defaultValue);
                const auto uTau = kappa * initialVelocity[0] / log(real_t(1) / roughnessLengthRatio);
                pressureGradient_ = uTau * uTau / real_t(domainSetup.domainSize_[2]);
                windDirection_ = parameters.getParameter<Value_T>("windDirectionDegrees", real_t(0.0));
                windDirection_ *= walberla::math::pi / real_t(180);

                // Dynamic wind-direction pressure-driven
                windDirection_ = parameters.getParameter<Value_T>("initialWindDirectionDegrees", real_t(0.0));

                const std::string windDirChangeRateTable = parameters.getParameter<std::string>("windDirectionChangeRateTable", "");

                if(!windDirChangeRateTable.empty()) {
                    std::ifstream file(windDirChangeRateTable);
                    std::string line;

                    std::vector<real_t> times;
                    std::vector<real_t> windDirChangeRate;

                    if (!file.is_open()) {
                        WALBERLA_ABORT("Error opening file " << windDirChangeRateTable << ", please check file existence.");
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
                            real_t time, changeRate;
                            in >> time >> changeRate;
                            times.push_back(time  / Conversion::C_t());                    // s as an input
                            windDirChangeRate.emplace_back(changeRate * Conversion::C_t());  // rad/s as an input
                        }

                        file.close();
                    }

                    windDirectionChangeInterpolator_.setPoints(times.size(), times.data(), windDirChangeRate.data());
                    timesTable_ = times;
                    windDirChangeRateTable_ = windDirChangeRate;
                }

                // Coriolis-driven
                coriolisFrequency_ = parameters.getParameter<Value_T>("coriolisFrequency", defaultValue);
                coriolisFrequency_ *= Conversion::C_t();
                geostrophicWind_ = parameters.getParameter<walberla::Vector3<Value_T>>("geostrophicWind", walberla::Vector3<Value_T>(defaultValue));
                geostrophicWind_ *= (Conversion::C_t() / Conversion::C_l());

            }

            void operator()(walberla::IBlock * block, const uint_t level = 0, const uint_t executionCount = 0)
            {
                if( !walberla::selectable::isSetSelected( block->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
                    return;

                const uint_t currentIteration = timeloop_->getCurrentTimeStep();

                static_cast<T*>(this)->driveFlow(block, level, executionCount);
            }

            void operator()(const uint_t level = 0, const uint_t executionCount = 0) {
                for( auto it = forest_->begin(); it != forest_->end(); ++it ) {
                    this->operator()(it.get(), level, executionCount);
                }
            }

        protected:

            const BlockDataID forceFieldId_{};
            const BlockDataID velocityFieldId_{};

            DriverSetup::Type driverType_{};
            std::shared_ptr<Driver> driver_{};

            walberla::timeloop::ITimeloop * timeloop_{nullptr};

            Value_T pressureGradient_{};
            Value_T windDirection_{};
            Value_T coriolisFrequency_{};
            data_interpolation::DataInterpolator<Value_T, Value_T, false> windDirectionChangeInterpolator_;

            std::vector<real_t> timesTable_;
            std::vector<real_t> windDirChangeRateTable_;

            walberla::Vector3<Value_T> geostrophicWind_{};

        private:

            virtual void driveFlow(walberla::IBlock * block, const uint_t level, const uint_t executionCount) const {
                if(driver_)
                    driver_->driveFlow(block, level, executionCount);
            }

            const std::shared_ptr<walberla::blockforest::StructuredBlockForest> forest_{};

            const walberla::Set<walberla::SUID> requiredBlockSelectors_;
            const walberla::Set<walberla::SUID> incompatibleBlockSelectors_;
        };

        template<typename VectorField_T>
        struct FlowDriver : public FlowDriverBase<FlowDriver<VectorField_T>, typename VectorField_T::value_type> {

            friend FlowDriverBase<FlowDriver<VectorField_T>, typename VectorField_T::value_type>;

            FlowDriver( const std::shared_ptr<walberla::blockforest::StructuredBlockForest> & forest,
                        std::shared_ptr<walberla::Config> & globalConfig,
                        const domain::DomainSetup & domainSetup,
                        const BlockDataID forceFieldID, const BlockDataID velocityFieldID,
                        walberla::timeloop::ITimeloop * timeloop,
                        const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                        const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
            : FlowDriverBase<FlowDriver<VectorField_T>, typename VectorField_T::value_type>(
                    forest, globalConfig, domainSetup, forceFieldID, velocityFieldID, timeloop, requiredBlockSelectors, incompatibleBlockSelectors)
            {
                switch (this->driverType_) {
                    case DriverSetup::None :
                        this->driver_ = std::make_shared<ZeroForce<VectorField_T>>(this->forceFieldId_); break;
                    case DriverSetup::PressureGradient :
                        this->driver_ = std::make_shared<PressureGradient<VectorField_T>>(this->forceFieldId_, this->pressureGradient_, this->windDirection_); break;
                    case DriverSetup::DynamicPressureGradient :
                        this->driver_ = std::make_shared<DynamicPressureGradient<VectorField_T>>(this->velocityFieldId_, this->forceFieldId_,
                            this->pressureGradient_, this->windDirection_, this->timeloop_, this->windDirectionChangeInterpolator_); break;
                    case DriverSetup::CoriolisForce :
                        this->driver_ = std::make_shared<CoriolisForce<VectorField_T>>(this->velocityFieldId_, this->forceFieldId_, this->coriolisFrequency_,
                            this->geostrophicWind_); break;
                    default:
                        WALBERLA_ABORT("Invalid driver type: " << this->driverType_)
                }
            }

        };

    } // namespace external_forces

} // namespace turbine_core

#endif // TURBINECORE_FLOWDRIVER_H