
#pragma once

#ifndef TURBINECORE_FORCEMODEL_H
#define TURBINECORE_FORCEMODEL_H

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

#include "walberla_helper/field/Field.h"
#include "walberla_helper/field/Projectors.h"

#include "point3/Point3.h"
#include "point3/ActuatorData.h"

#include "walberla_helper/blockforest/BlockInfo.h"

namespace turbine_core {

    namespace force_model {

        class ForceModel {

        public:

            HOST_DEVICE_PREFIX ForceModel() {}

            HOST_DEVICE_PREFIX virtual ~ForceModel() {}

            HOST_DEVICE_PREFIX void print() const {
                do_print();
            }

            HOST_DEVICE_PREFIX ForceModel * clone() const {
                return do_clone();
            }

            HOST_DEVICE_PREFIX void evaluateDensityAndVelocity( const blockforest::BlockInfo & blockInfo,
                                                                const projectors::Projectors & projectors,
                                                                const Point3 <real_t> *points, uint_t pointIdy) {
                do_evaluateDensityAndVelocity(blockInfo, projectors, points, pointIdy);
            }

            HOST_DEVICE_PREFIX void calculateForces( const Point3<real_t> * points,
                                                     const Point3<real_t> * parentPoint) {
                do_calculateForces(points, parentPoint);
            }

            HOST_DEVICE_PREFIX void spreadForces( const blockforest::BlockInfo & blockInfo,
                                                  const projectors::Projectors & projectors,
                                                  const Point3 <real_t> *points, uint_t pointIdy) {
                do_spreadForces(blockInfo, projectors, points, pointIdy);
            }

            HOST_DEVICE_PREFIX void resetCommunicationData() { do_resetCommunicationData(); }

            HOST_DEVICE_PREFIX virtual ActuatorData * points() {
                printf("Base class ForceModel is not implemented.\n");
                exit(EXIT_FAILURE);
            }

            HOST_DEVICE_PREFIX virtual bool * remotePoint() {
                printf("Base class ForceModel is not implemented.\n");
                exit(EXIT_FAILURE);
            }

        private:

            HOST_DEVICE_PREFIX virtual void do_print() const {
                printf("Base class ForceModel is not implemented.\n");
                exit(EXIT_FAILURE);
            }

            HOST_DEVICE_PREFIX virtual ForceModel * do_clone() const {
                printf("Base class ForceModel is not implemented.\n");
                exit(EXIT_FAILURE);
            }

            HOST_DEVICE_PREFIX virtual void do_evaluateDensityAndVelocity( const blockforest::BlockInfo & blockInfo,
                                                                           const projectors::Projectors & projectors,
                                                                           const Point3 <real_t> *points,
                                                                           uint_t pointIdy) {
                exit(EXIT_FAILURE);
            }

            HOST_DEVICE_PREFIX virtual void do_calculateForces( const Point3<real_t> * points,
                                                                const Point3<real_t> * parentPoint) {
                exit(EXIT_FAILURE);
            }

            HOST_DEVICE_PREFIX virtual void do_spreadForces( const blockforest::BlockInfo & blockInfo,
                                                             const projectors::Projectors & projectors,
                                                             const Point3 <real_t> *points,
                                                             uint_t pointIdy) {
                exit(EXIT_FAILURE);
            }

            HOST_DEVICE_PREFIX virtual void do_resetCommunicationData() {}
        };
    }

    using force_model::ForceModel;

}

#endif //TURBINECORE_FORCEMODEL_H
