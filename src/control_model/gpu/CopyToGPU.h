
#pragma once

#ifndef TURBINECORE_CONTROLMODEL_COPYTOGPU_H
#define TURBINECORE_CONTROLMODEL_COPYTOGPU_H

#include <cassert>
#include <vector>

#include "wind_turbine_core/ProjectDefines.h"

#include "point3/ActuatorData.h"

#include "control_model/all.h"

namespace turbine_core {

    namespace control_model {

        namespace gpu {

            GLOBAL_PREFIX void deleteControlModelOnGPU( ControlModel ** gpuPtr ) {
                delete (*gpuPtr);
            }

            template<class ControlModel_T>
            struct CopyToGPU {

                HOST_PREFIX CopyToGPU() {

                    cudaMalloc((void**)&devicePtr, sizeof(ControlModel_T*));
                    gpuErrchk(cudaPeekAtLastError())
                    gpuErrchk(cudaDeviceSynchronize())

                }

                HOST_PREFIX ~CopyToGPU() {
                    deleteControlModelOnGPU<<<1,1>>>((ControlModel**)devicePtr);
                    gpuErrchk(cudaPeekAtLastError())
                    gpuErrchk(cudaDeviceSynchronize())
                    cudaFree(devicePtr);
                }

                ControlModel_T ** devicePtr{nullptr};

            };

            // template prototypes
            template<typename ControlModel_T>
            HOST_PREFIX void copy( CopyToGPU<ControlModel_T> &, const std::shared_ptr<ControlModel_T> & );

            GLOBAL_PREFIX void newRAWSAnglesModelOnGPU( RAWSAnglesModel ** gpuPtr,
                const real_t alphaFilter, const real_t previousRAWS, const real_t angleX,
                    const real_t angleY, const real_t angleZ, const uint_t nPoints,
                    real_t * windVelocities, Vector3<real_t> * angles)
            {
                (*gpuPtr) = new RAWSAnglesModel(alphaFilter, previousRAWS, angleX, angleY,
                                                angleZ, nPoints, windVelocities, angles);
            }

            GLOBAL_PREFIX void newRAWSOmegaModelOnGPU( RAWSOmegaModel ** gpuPtr, const real_t alphaFilter,
                                                       const real_t previousRAWS, const uint_t nPoints,
                                                       real_t * windVelocities, Vector3<real_t> * rotationalVelocities )
            {
                (*gpuPtr) = new RAWSOmegaModel(alphaFilter, previousRAWS, nPoints, windVelocities, rotationalVelocities);
            }

            GLOBAL_PREFIX void newTimeControlVelocitiesModel( TimeControlVelocitiesModel ** gpuPtr,
                                                              const uint_t nPoints,
                                                              real_t * times, Vector3<real_t> * translationalVelocities, Vector3<real_t> * rotationalVelocities )
            {
                (*gpuPtr) = new TimeControlVelocitiesModel(nPoints, times, translationalVelocities, rotationalVelocities);
            }

            // template specialisation
            GLOBAL_PREFIX void newDISCONTorqueModel( DISCONTorqueModel ** gpuPtr, const real_t cutinSpeed,
                    const real_t region_2_startingspeed, const real_t region_2_endingspeed, const real_t rated_generator_speed,
                    const real_t cutinTorque, const real_t ratedGeneratorTorque, const real_t rotorRadius, const real_t Ngear,
                    const real_t ratedPower, const real_t CpOpt, const real_t tsrOpt, const real_t iRot, const real_t iGen,
                    const real_t gearboxEfficiency, const real_t generatorEfficiency)
            {
                (*gpuPtr) = new DISCONTorqueModel(cutinSpeed, region_2_startingspeed, region_2_endingspeed, rated_generator_speed,
                    cutinTorque, ratedGeneratorTorque, rotorRadius, Ngear, ratedPower, CpOpt, tsrOpt, iRot, iGen, gearboxEfficiency, generatorEfficiency);
            }

            // template specialisation
            GLOBAL_PREFIX void newDISCONPitchModel( DISCONPitchModel ** gpuPtr, const real_t previousPitch,
                    const real_t pitchK, const real_t pitchControlKP, const real_t pitchControlKI, const real_t pitchMin,
                    const real_t pitchMax, const real_t ratedGeneratorSpeed, const real_t nGear, const real_t maxPitchRate)
            {
                (*gpuPtr) = new DISCONPitchModel(previousPitch, pitchK, pitchControlKP, pitchControlKI, pitchMin,
                    pitchMax, ratedGeneratorSpeed, nGear, maxPitchRate);
            }

            // template specialisation
            GLOBAL_PREFIX void newSimpleYawController( SimpleYawController ** gpuPtr, const real_t yawingSpeed)
            {
                (*gpuPtr) = new SimpleYawController(yawingSpeed);
            }

            template<> HOST_PREFIX void copy( CopyToGPU<ControlModel>& cpyToGPU, const std::shared_ptr<ControlModel> & cpuPtr ) {
                WALBERLA_ABORT("ControlModel: abstract base class cannot be copied.")
            }

            // template specialisation
            template<> HOST_PREFIX void copy( CopyToGPU<RAWSAnglesModel>& cpyToGPU, const std::shared_ptr<RAWSAnglesModel> & cpuPtr ) {

                // windVelocity and omegas are a vector of real_t and a vector of vector3<real_t>
                // we get the size of the vector
                uint_t nPoints = cpuPtr->interpolator()->length();

                // create pointers to device data and allocate
                real_t *d_windVelocity;
                Vector3<real_t> *d_angles;
                cudaMalloc((void**)&d_windVelocity,  nPoints * sizeof(real_t));
                cudaMalloc((void**)&d_angles,        nPoints * sizeof(Vector3<real_t>));

                // copy the data from the cpu to the gpu
                cudaMemcpy(d_windVelocity,  cpuPtr->interpolator()->x(),  nPoints * sizeof(real_t), cudaMemcpyHostToDevice);
                cudaMemcpy(d_angles,  cpuPtr->interpolator()->y(),  nPoints * sizeof(Vector3<real_t>), cudaMemcpyHostToDevice);

                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

                newRAWSAnglesModelOnGPU<<<1,1>>>(cpyToGPU.devicePtr, cpuPtr->alphaFilter(), cpuPtr->previousRAWS(),
                                                cpuPtr->previousAngles()[0], cpuPtr->previousAngles()[1], cpuPtr->previousAngles()[2],
                                                nPoints, d_windVelocity, d_angles);

                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

                cudaFree(d_windVelocity);
                cudaFree(d_angles);

                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

            }

            // template specialisation
            template<> HOST_PREFIX void copy( CopyToGPU<RAWSOmegaModel>& cpyToGPU, const std::shared_ptr<RAWSOmegaModel> & cpuPtr ) {

                // windVelocity and omegas are a vector of real_t and a vector of vector3<real_t>
                // we get the size of the vector
                uint_t nPoints = cpuPtr->interpolator()->length();

                // create pointers to device data and allocate
                real_t *d_windVelocity;
                Vector3<real_t> *d_omegas;
                cudaMalloc((void**)&d_windVelocity,  nPoints * sizeof(real_t));
                cudaMalloc((void**)&d_omegas,        nPoints * sizeof(Vector3<real_t>));

                // copy the data from the cpu to the gpu
                cudaMemcpy(d_windVelocity,  cpuPtr->interpolator()->x(),  nPoints * sizeof(real_t), cudaMemcpyHostToDevice);
                cudaMemcpy(d_omegas,  cpuPtr->interpolator()->y(),  nPoints * sizeof(Vector3<real_t>), cudaMemcpyHostToDevice);

                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

                newRAWSOmegaModelOnGPU<<<1,1>>>(cpyToGPU.devicePtr, cpuPtr->alphaFilter(), cpuPtr->previousRAWS(), nPoints, d_windVelocity, d_omegas);

                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

                cudaFree(d_windVelocity);
                cudaFree(d_omegas);

                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

            }

            // template specialisation
            template<> HOST_PREFIX void copy( CopyToGPU<TimeControlVelocitiesModel>& cpyToGPU, const std::shared_ptr<TimeControlVelocitiesModel> & cpuPtr ) {

                // windVelocity and omegas are a vector of real_t and a vector of vector3<real_t>
                // we get the size of the vector
                uint_t nPoints = cpuPtr->interpolatorTranslations()->length();

                // create pointers to device data and allocate
                real_t *d_times;
                Vector3<real_t> *d_translation;
                Vector3<real_t> *d_rotations;
                cudaMalloc((void**)&d_times,  nPoints * sizeof(real_t));
                cudaMalloc((void**)&d_translation,  nPoints * sizeof(Vector3<real_t>));
                cudaMalloc((void**)&d_rotations,    nPoints * sizeof(Vector3<real_t>));

                // copy the data from the cpu to the gpu
                cudaMemcpy(d_times,  cpuPtr->interpolatorTranslations()->x(),  nPoints * sizeof(real_t), cudaMemcpyHostToDevice);
                cudaMemcpy(d_translation,  cpuPtr->interpolatorTranslations()->y(),  nPoints * sizeof(Vector3<real_t>), cudaMemcpyHostToDevice);
                cudaMemcpy(d_rotations,  cpuPtr->interpolatorRotations()->y(),  nPoints * sizeof(Vector3<real_t>), cudaMemcpyHostToDevice);

                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

                newTimeControlVelocitiesModel<<<1,1>>>(cpyToGPU.devicePtr, nPoints, d_times, d_translation, d_rotations);

                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

                cudaFree(d_times);
                cudaFree(d_translation);
                cudaFree(d_rotations);

                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

            }

            // template specialisation
            template<> HOST_PREFIX void copy( CopyToGPU<DISCONTorqueModel>& cpyToGPU, const std::shared_ptr<DISCONTorqueModel> & cpuPtr ) {

                newDISCONTorqueModel<<<1,1>>>(cpyToGPU.devicePtr, cpuPtr->cutinSpeed(), cpuPtr->region_2_startingspeed(), cpuPtr->region_2_endingspeed(),
                    cpuPtr->rated_generator_speed(), cpuPtr->cutinTorque(), cpuPtr->ratedGeneratorTorque(), cpuPtr->rotorRadius(), cpuPtr->Ngear(),
                    cpuPtr->ratedPower(), cpuPtr->CpOpt(), cpuPtr->tsrOpt(), cpuPtr->iRot(), cpuPtr->iGen(), cpuPtr->gearboxEfficiency(),
                    cpuPtr->generatorEfficiency()
                );

                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())
            }

            // template specialisation
            template<> HOST_PREFIX void copy( CopyToGPU<DISCONPitchModel>& cpyToGPU, const std::shared_ptr<DISCONPitchModel> & cpuPtr ) {

                newDISCONPitchModel<<<1,1>>>(cpyToGPU.devicePtr, cpuPtr->previousPitch(), cpuPtr->pitchK(), cpuPtr->pitchControlKP(), cpuPtr->pitchControlKI(),
                    cpuPtr->pitchMin(), cpuPtr->pitchMax(), cpuPtr->ratedGeneratorSpeed(), cpuPtr->nGear(), cpuPtr->maxPitchRate()
                );

                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())
            }

            // template specialisation
            template<> HOST_PREFIX void copy( CopyToGPU<SimpleYawController>& cpyToGPU, const std::shared_ptr<SimpleYawController> & cpuPtr ) {

                newSimpleYawController<<<1,1>>>(cpyToGPU.devicePtr, cpuPtr->yawingSpeed());

                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())
            }
        }

    }

}

#endif //TURBINECORE_CONTROLMODEL_COPYTOGPU_H
