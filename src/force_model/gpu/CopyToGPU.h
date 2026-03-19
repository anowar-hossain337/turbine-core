
#pragma once

#ifndef TURBINECORE_FORCEMODEL_COPYTOGPU_H
#define TURBINECORE_FORCEMODEL_COPYTOGPU_H

#include <cassert>
#include <vector>

#include "wind_turbine_core/ProjectDefines.h"

#include "point3/ActuatorData.h"

namespace turbine_core {

    namespace force_model {

        class ForceModel;
        class None;
        template< typename DensityAndVelocityInterpolator_T, typename ForceDistributor_T, typename TipLossModel_T >
        class ActuatorLineModel;
        template< typename DensityAndVelocityInterpolator_T, typename ForceDistributor_T >
        class ActuatorDiskModel;

        class TipLossModel;

        namespace gpu {

            GLOBAL_PREFIX void deleteForceModelOnGPU( ForceModel ** gpuPtr ) {
                delete (*gpuPtr);
            }

            template<class ForceModel_T>
            struct CopyToGPU {

                HOST_PREFIX CopyToGPU() {

                    cudaMalloc((void**)&devicePtr, sizeof(ForceModel_T*));
                    gpuErrchk(cudaPeekAtLastError())
                    gpuErrchk(cudaDeviceSynchronize())

                }

                HOST_PREFIX ~CopyToGPU() {
                    deleteForceModelOnGPU<<<1,1>>>((ForceModel**)devicePtr);
                    gpuErrchk(cudaPeekAtLastError())
                    gpuErrchk(cudaDeviceSynchronize())
                    cudaFree(devicePtr);
                }

                ForceModel_T ** devicePtr{nullptr};

            };

            // template prototypes
            template<typename ForceModel_T>
            HOST_PREFIX void copy( CopyToGPU<ForceModel_T> &, const std::shared_ptr<ForceModel_T> & );

            template<typename DensityAndVelocityInterpolator_T, typename ForceDistributor_T, template<typename...> class ForceModel_T>
            HOST_PREFIX void copy( CopyToGPU<ForceModel_T<DensityAndVelocityInterpolator_T,ForceDistributor_T>> &,
                                   const std::shared_ptr<ForceModel_T<DensityAndVelocityInterpolator_T,ForceDistributor_T>> & );


            // kernel

            template< typename DensityAndVelocityInterpolator_T, typename ForceDistributor_T, typename TipLossModel_T >
            GLOBAL_PREFIX void newActuatorLineModelOnGPU( ActuatorLineModel<DensityAndVelocityInterpolator_T,ForceDistributor_T,TipLossModel_T> ** gpuPtr, const real_t elementLength,
                                                          const uint_t nPoints, ActuatorData * _data, uint_t * nInterpolationPoints, real_t ** aoa, Vector2<real_t> ** polar,
                                                          TipLossModel_T* tipLossModel)
            {
                //TODO might want to think of something cleverer!
                auto data = (ActuatorData*) malloc(nPoints * sizeof(ActuatorData));
                for(uint_t i = 0; i < nPoints; ++i) {
                    new (static_cast<ActuatorData*>(data) + i) ActuatorData(_data[i].width, aerodynamics::AirfoilPolar(nInterpolationPoints[i], aoa[i], polar[i]));
                }

                (*gpuPtr) = new ActuatorLineModel<DensityAndVelocityInterpolator_T, ForceDistributor_T, TipLossModel_T>(elementLength, nPoints, data, *tipLossModel);

                for(uint_t i = 0; i < nPoints; ++i) {
                    (static_cast<ActuatorData*>(data) + i)->~ActuatorData();
                }
                free(data);
            }

            template< typename DensityAndVelocityInterpolator_T, typename ForceDistributor_T >
            GLOBAL_PREFIX void newActuatorDiskModelOnGPU( ActuatorDiskModel<DensityAndVelocityInterpolator_T,ForceDistributor_T> ** gpuPtr, const real_t radius, const real_t dragCoefficient )
            {
                (*gpuPtr) = new ActuatorDiskModel<DensityAndVelocityInterpolator_T,ForceDistributor_T>(radius, dragCoefficient);
            }

            // template specialisation
            template<> HOST_PREFIX void copy( CopyToGPU<ForceModel>& cpyToGPU, const std::shared_ptr<ForceModel> & cpuPtr ) {
                WALBERLA_ABORT("Abstract base class cannot be copied.")
            }

            // template specialisations

            template< typename DensityAndVelocityInterpolator_T, typename ForceDistributor_T, typename TipLossModel_T >
            HOST_PREFIX void copy( CopyToGPU<ActuatorLineModel<DensityAndVelocityInterpolator_T,ForceDistributor_T,TipLossModel_T>>& cpyToGPU,
                                   const std::shared_ptr<ActuatorLineModel<DensityAndVelocityInterpolator_T,ForceDistributor_T,TipLossModel_T>> & cpuPtr ) {

                uint_t nPoints = cpuPtr->nPoints_;

                // copy actuator data
                ActuatorData * d_actuatorData;
                cudaMalloc((void**)&d_actuatorData, nPoints * sizeof(ActuatorData) );
                cudaMemcpy(d_actuatorData, cpuPtr->points_, nPoints * sizeof(ActuatorData), cudaMemcpyHostToDevice);

                // copy airfoil polar of actuator data
                // must be handled separately as it contains dynamic memory
                uint_t * d_nInterpolationPoints;
                cudaMalloc((void**)&d_nInterpolationPoints, nPoints * sizeof(uint_t));

                real_t **d_aoa;
                Vector2<real_t> **d_polar;
                cudaMalloc((void**)&d_aoa,  nPoints * sizeof(real_t*));
                cudaMalloc((void**)&d_polar,  nPoints * sizeof(Vector2<real_t>*));

                std::vector<uint_t> nInterpolationPoints(nPoints);

                auto **h_d_aoa = (real_t**)malloc(nPoints * sizeof(real_t*));
                auto **h_d_polar = (Vector2<real_t>**)malloc(nPoints * sizeof(Vector2<real_t>*));

                for(uint_t i = 0; i < nPoints; ++i) {

                    nInterpolationPoints[i] = cpuPtr->points_[i].polar.length();

                    cudaMalloc((void **) &h_d_aoa[i],  nInterpolationPoints[i] * sizeof(real_t));
                    cudaMalloc((void **) &h_d_polar[i],  nInterpolationPoints[i] * sizeof(Vector2<real_t>));

                    cudaMemcpy(h_d_aoa[i], cpuPtr->points_[i].polar.interpolator.x(), nInterpolationPoints[i] * sizeof(real_t), cudaMemcpyHostToDevice);
                    gpuErrchk(cudaPeekAtLastError())
                    gpuErrchk(cudaDeviceSynchronize())

                    cudaMemcpy(h_d_polar[i], cpuPtr->points_[i].polar.interpolator.y(), nInterpolationPoints[i] * sizeof(Vector2<real_t>), cudaMemcpyHostToDevice);
                    gpuErrchk(cudaPeekAtLastError())
                    gpuErrchk(cudaDeviceSynchronize())

                }

                cudaMemcpy(d_aoa,  h_d_aoa,  nPoints * sizeof(real_t*), cudaMemcpyHostToDevice);
                cudaMemcpy(d_polar,  h_d_polar,  nPoints * sizeof(Vector2<real_t>*), cudaMemcpyHostToDevice);
                cudaMemcpy(d_nInterpolationPoints, nInterpolationPoints.data(), nPoints * sizeof(uint_t), cudaMemcpyHostToDevice);

                // copy tip loss model
                TipLossModel_T * d_tipLossModel;
                cudaMalloc((void**)&d_tipLossModel, sizeof(TipLossModel_T) );
                cudaMemcpy(d_tipLossModel, &cpuPtr->tiplossModel_, sizeof(TipLossModel_T), cudaMemcpyHostToDevice);

                newActuatorLineModelOnGPU<<<1,1>>>(cpyToGPU.devicePtr, cpuPtr->elementLength(), nPoints, d_actuatorData, d_nInterpolationPoints, d_aoa, d_polar, d_tipLossModel);
                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

                cudaFree(d_actuatorData);
                cudaFree(d_nInterpolationPoints);

                for(uint_t i = 0; i < nPoints; ++i) {
                    cudaFree(h_d_aoa[i]);
                    cudaFree(h_d_polar[i]);
                }

                free(h_d_aoa);
                free(h_d_polar);

                cudaFree(d_aoa);
                cudaFree(d_polar);

                cudaFree(d_tipLossModel);

                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())
            }

            template< typename DensityAndVelocityInterpolator_T, typename ForceDistributor_T >
            HOST_PREFIX void copy( CopyToGPU<ActuatorDiskModel<DensityAndVelocityInterpolator_T,ForceDistributor_T>>& cpyToGPU,
                                   const std::shared_ptr<ActuatorDiskModel<DensityAndVelocityInterpolator_T,ForceDistributor_T>> & cpuPtr ) {

                uint_t nPoints = 1;

                // copy actuator data
                // ActuatorData * d_actuatorData;
                // cudaMalloc((void**)&d_actuatorData, nPoints * sizeof(ActuatorData) );
                // cudaMemcpy(d_actuatorData, cpuPtr->point_, nPoints * sizeof(ActuatorData), cudaMemcpyHostToDevice);

                newActuatorDiskModelOnGPU<<<1,1>>>(cpyToGPU.devicePtr, cpuPtr->radius_, cpuPtr->drag_);
                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

                // cudaFree(d_actuatorData);

            }
        }

    }

}

#endif //TURBINECORE_FORCEMODEL_COPYTOGPU_H
