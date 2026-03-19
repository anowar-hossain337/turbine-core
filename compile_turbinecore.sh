#!/bin/bash

cmake .. -DWALBERLA_DIR=/home/yanlab/walberla/walberla-dev \
         -DPython_ROOT_DIR=/home/yanlab/walberla/venv-walberla/bin \
         -DWALBERLA_BUILD_WITH_CODEGEN=ON \
         -DWALBERLA_BUILD_WITH_PYTHON=ON
