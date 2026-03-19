cmake_minimum_required(VERSION 3.15)

project(walberla-wind
        VERSION 1.0.0
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find waLBerla
find_package(waLBerla REQUIRED CONFIG)
include(${WALBERLA_DIR}/CMake/walberla_macros.cmake)

# Optional MPI support (inherited from waLBerla)
if(WALBERLA_BUILD_WITH_MPI)
    find_package(MPI REQUIRED)
endif()

enable_testing()

add_subdirectory(src)
add_subdirectory(apps)
add_subdirectory(tests)
