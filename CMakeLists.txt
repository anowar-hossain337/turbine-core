CMAKE_MINIMUM_REQUIRED (VERSION 3.24)

if(WALBERLA_BUILD_WITH_CUDA)
    PROJECT ( TurbineCore LANGUAGES CXX CUDA )
    if (CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 12.0)
        message(FATAL_ERROR "CUDA version must be at least 12.0!  " ${CMAKE_CUDA_COMPILER_VERSION})
    endif()

    find_package(CUDAToolkit REQUIRED)
    set( CMAKE_CUDA_SEPARABLE_COMPILATION ON )

    if (CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 11.0)
        set(CMAKE_CUDA_STANDARD 14)
        set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)
    endif()

else()
    PROJECT ( TurbineCore LANGUAGES CXX )
endif()

#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -rdc=true -c")

# Extends cmake module path - so that FindwaLBerla.cmake in the current directory is found
set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${TurbineCore_SOURCE_DIR} )
find_package( waLBerla )

enable_testing()
include_directories( src )
include_directories ( ${wind_turbine_BINARY_DIR}/src )

# enable waLBerla-wind Python modules
set( WALBERLA_WIND_PYTHON_DIR ${TurbineCore_SOURCE_DIR}/python )
set( ENV{PYTHONPATH} "${WALBERLA_WIND_PYTHON_DIR}:$ENV{PYTHONPATH}" )

# Simple compatibility wrapper: some apps expect a helper
# function `waLBerla_add_executable(NAME <name> FILES <files> DEPENDS <libs...>)`.
# waLBerla (upstream) provides various helper functions; if it's missing
# here we provide a minimal implementation that maps to CMake's add_executable
# and target_link_libraries. This avoids the "Unknown CMake command"
# error when configuring the apps folder.
function(waLBerla_add_executable)
    set(oneValueArgs NAME)
    set(multiValueArgs FILES DEPENDS)
    cmake_parse_arguments(WALEXE "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    if(NOT WALEXE_NAME)
        message(FATAL_ERROR "waLBerla_add_executable called without NAME")
    endif()

    # create the target
    add_executable(${WALEXE_NAME} ${WALEXE_FILES})

    # link requested dependencies (if any)
    if(WALEXE_DEPENDS)
        target_link_libraries(${WALEXE_NAME} PUBLIC ${WALEXE_DEPENDS})
    endif()

    # make the target available to parent scopes
    set_target_properties(${WALEXE_NAME} PROPERTIES FOLDER "apps")
endfunction()

# Minimal compatibility for waLBerla_compile_test used by turbinecore tests.
# It creates an executable target for each provided source file. The target
# name is derived from the filename without extension. This is intentionally
# small: it allows CMake configuration to succeed. Linking against specific
# walberla libraries can be added later if needed for building the tests.
function(waLBerla_compile_test)
    set(oneValueArgs)
    set(multiValueArgs FILES)
    cmake_parse_arguments(WACTEST "" "" "${multiValueArgs}" ${ARGN})

    if(NOT WACTEST_FILES)
        message(FATAL_ERROR "waLBerla_compile_test called without FILES")
    endif()

    foreach(src ${WACTEST_FILES})
        # derive a target name from file name (strip directory and extension)
        get_filename_component(fname ${src} NAME_WE)
        # create the executable target in the current directory
        add_executable(${fname} ${src})
        # make it visible in IDEs under tests
        set_target_properties(${fname} PROPERTIES FOLDER "tests")
    endforeach()
endfunction()

add_subdirectory( src  )
add_subdirectory( apps )

if ( WALBERLA_BUILD_TESTS )
    add_subdirectory ( tests )
else()
    add_subdirectory( tests EXCLUDE_FROM_ALL )
endif()
