#pragma once

#ifndef TURBINECORE_DISCRETISATION_COPYTOGPU_H
#define TURBINECORE_DISCRETISATION_COPYTOGPU_H

#include <cassert>

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/WalberlaDataTypes.h"

#include "point3/Point3.h"

#include "discretisation/all.h"

namespace turbine_core {

    namespace discretisation {

        namespace gpu {

            GLOBAL_PREFIX void deleteDiscretisationOnGPU( Discretisation ** gpuPtr ) {
                delete (*gpuPtr);
            }

            template<typename Discretisation_T>
            struct CopyToGPU {

                HOST_PREFIX CopyToGPU() {

                    cudaMalloc((void**)&devicePtr, sizeof(Discretisation_T*));
                    gpuErrchk(cudaPeekAtLastError())
                    gpuErrchk(cudaDeviceSynchronize())

                }

                HOST_PREFIX ~CopyToGPU() {
                    deleteDiscretisationOnGPU<<<1,1>>>((Discretisation**)devicePtr);
                    gpuErrchk(cudaPeekAtLastError())
                    gpuErrchk(cudaDeviceSynchronize())
                    cudaFree(devicePtr);
                }

                Discretisation_T ** devicePtr{nullptr};

            };

            // kernel
            GLOBAL_PREFIX void newDiskDiscretisationOnGPU( Disk ** gpuPtr, const Disk * disk ) {
                (*gpuPtr) = new Disk(*disk);
            }


            GLOBAL_PREFIX void newLineDiscretisationOnGPU( Line ** gpuPtr, const Line * line, const uint_t nPoints, Point3<real_t> * points ) {
                (*gpuPtr) = new Line(*line, nPoints, points);
            }


            // template prototype
            template<typename Discretisation_T>
            HOST_PREFIX void copy( CopyToGPU<Discretisation_T> &, const std::shared_ptr<Discretisation_T> & );


            // template specialisation
            template<> HOST_PREFIX void copy( CopyToGPU<Disk>& cpyToGPU, const std::shared_ptr<Disk> & cpuPtr ) {

                Disk * d_disk;
                cudaMalloc((void**)&d_disk, sizeof(Disk));
                cudaMemcpy(d_disk, cpuPtr.get(), sizeof(Disk), cudaMemcpyHostToDevice);
                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

                newDiskDiscretisationOnGPU<<<1,1>>>(cpyToGPU.devicePtr, d_disk);
                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

                cudaFree(d_disk);

            }

            template<> HOST_PREFIX void copy( CopyToGPU<Line>& cpyToGPU, const std::shared_ptr<Line> & cpuPtr ) {

                Line * d_line;
                cudaMalloc((void**)&d_line, sizeof(Line));
                cudaMemcpy(d_line, cpuPtr.get(), sizeof(Line), cudaMemcpyHostToDevice);
                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

                Point3<real_t> * d_points;
                cudaMalloc((void**)&d_points, cpuPtr->nPoints() * sizeof(Point3<real_t>));
                cudaMemcpy(d_points, cpuPtr->points(), cpuPtr->nPoints() * sizeof(Point3<real_t>), cudaMemcpyHostToDevice);
                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

                newLineDiscretisationOnGPU<<<1,1>>>(cpyToGPU.devicePtr, d_line, cpuPtr->nPoints(), d_points);
                gpuErrchk(cudaPeekAtLastError())
                gpuErrchk(cudaDeviceSynchronize())

                cudaFree(d_line);
                cudaFree(d_points);
            }


        } // namespace gpu

    } // namespace discretisation

} // namespace turbine_core

#endif //TURBINECORE_DISCRETISATION_COPYTOGPU_H
