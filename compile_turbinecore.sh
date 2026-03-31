#!/bin/bash

cmake .. -DWALBERLA_DIR=/home/yanlab/walberla/walberla-dev \
         -DPython_ROOT_DIR=/home/yanlab/walberla/venv-walberla/bin \
         -DWALBERLA_BUILD_WITH_CODEGEN=ON \
         -DWALBERLA_BUILD_WITH_PYTHON=ON
         
cmake \
  -DWALBERLA_DIR=/home/yanpersonalresearch/walberla/walberla-wind \
  -DPython_ROOT_DIR=/home/yanpersonalresearch/walberla/myenv \
  -DWALBERLA_BUILD_WITH_CODEGEN=ON \
  -DWALBERLA_BUILD_WITH_PYTHON=ON \
  -DWALBERLA_BUILD_WITH_CUDA=ON \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda-13.2/bin/nvcc \
  -DCUDAToolkit_ROOT=/usr/local/cuda-13.2 \
  ..
  
  Error msg broken g++
  cmake \
  -DCMAKE_C_COMPILER=/usr/bin/gcc \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda-13.2/bin/nvcc \
  -DCUDAToolkit_ROOT=/usr/local/cuda-13.2 \
  -DWALBERLA_DIR=/home/yanpersonalresearch/walberla/walberla-wind \
  -DPython_ROOT_DIR=/home/yanpersonalresearch/walberla/myenv \
  -DWALBERLA_BUILD_WITH_CODEGEN=ON \
  -DWALBERLA_BUILD_WITH_PYTHON=ON \
  -DWALBERLA_BUILD_WITH_CUDA=ON \
  ../turbine-core
  
  error locating walberla
  cmake \
  -DCMAKE_C_COMPILER=/usr/bin/gcc \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda-13.2/bin/nvcc \
  -DCUDAToolkit_ROOT=/usr/local/cuda-13.2 \
  -DWALBERLA_DIR=/home/yanpersonalresearch/walberla/walberla \
  -DPython_ROOT_DIR=/home/yanpersonalresearch/walberla/myenv \
  -DWALBERLA_BUILD_WITH_CODEGEN=ON \
  -DWALBERLA_BUILD_WITH_PYTHON=ON \
  -DWALBERLA_BUILD_WITH_CUDA=ON \
  ../turbine-core
  
  New virtual environment has been built
  
    cmake \
  -DCMAKE_C_COMPILER=/usr/bin/gcc \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda-13.2/bin/nvcc \
  -DCUDAToolkit_ROOT=/usr/local/cuda-13.2 \
  -DWALBERLA_DIR=/home/yanpersonalresearch/walberla/walberla \
  -DPython_ROOT_DIR=/home/yanpersonalresearch/walberla/myenv \
  -DWALBERLA_BUILD_WITH_CODEGEN=ON \
  -DWALBERLA_BUILD_WITH_PYTHON=ON \
  -DWALBERLA_BUILD_WITH_CUDA=ON \
  ../turbine-core
