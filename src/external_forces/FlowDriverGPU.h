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
//! \file FlowDriverGPU.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef TURBINECORE_FLOWDRIVERGPU_H
#define TURBINECORE_FLOWDRIVERGPU_H

#include <functional>

#include <core/Set.h>
#include <core/selectable/IsSetSelected.h>

#include <core/uid/SUID.h>

#include <domain_decomposition/IBlock.h>
#include <blockforest/BlockForest.h>
#include <blockforest/StructuredBlockForest.h>

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "external_forces/FlowDriver.h"
#include "external_forces/DriverSetup.h"
#include "external_forces/ZeroForceGPU.h"
#include "external_forces/PressureGradientGPU.h"
#include "external_forces/DynamicPressureGradientGPU.h"
#include "external_forces/CoriolisForceGPU.h"

namespace turbine_core {

    namespace external_forces {

        template<typename GPUField_T>
        struct FlowDriverGPU : public FlowDriverBase<FlowDriverGPU<GPUField_T>, typename GPUField_T::value_type> {

            friend FlowDriverBase<FlowDriverGPU<GPUField_T>, typename GPUField_T::value_type>;

            FlowDriverGPU( const std::shared_ptr<walberla::blockforest::StructuredBlockForest> & forest,
                           std::shared_ptr<walberla::Config> & globalConfig,
                           const domain::DomainSetup & domainSetup,
                           const BlockDataID forceFieldID, const BlockDataID velocityFieldID,
                           walberla::timeloop::ITimeloop * timeloop,
                           const walberla::Set<walberla::SUID> & requiredBlockSelectors = walberla::Set<walberla::SUID>::emptySet(),
                           const walberla::Set<walberla::SUID> & incompatibleBlockSelectors = walberla::Set<walberla::SUID>::emptySet() )
            : FlowDriverBase<FlowDriverGPU<GPUField_T>, typename GPUField_T::value_type>(
                    forest, globalConfig, domainSetup, forceFieldID, velocityFieldID, timeloop, requiredBlockSelectors, incompatibleBlockSelectors)
            {
                switch (this->driverType_) {
                    case DriverSetup::None :
                        this->driver_ = std::make_shared<ZeroForceGPU<GPUField_T>>(this->forceFieldId_); break;
                    case DriverSetup::PressureGradient :
                        this->driver_ = std::make_shared<PressureGradientGPU<GPUField_T>>(this->forceFieldId_, this->pressureGradient_,
                        this->windDirection_); break;
                    case DriverSetup::DynamicPressureGradient :
                        this->driver_ = std::make_shared<DynamicPressureGradientGPU<GPUField_T>>(this->velocityFieldId_, this->forceFieldId_,
                            this->pressureGradient_, this->windDirection_, this->timeloop_, this->windDirectionChangeInterpolator_); break;
                    case DriverSetup::CoriolisForce :
                        this->driver_ = std::make_shared<CoriolisForceGPU<GPUField_T>>(this->velocityFieldId_, this->forceFieldId_, this->coriolisFrequency_, this->geostrophicWind_); break;
                    default:
                        WALBERLA_ABORT("Invalid driver type: " << this->driverType_)
                }
            }

        };

    } // namespace external_forces

} // namespace turbine_core

#endif // TURBINECORE_FLOWDRIVERGPU_H