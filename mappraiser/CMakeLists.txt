# MAPPRAISER

# Top-level CMakeLists file
# - contains the call to the project() command
# - set language standards and general options (e.g. PIC)
# - search for external dependencies: MPI, MKL, BLAS/LAPACK, ...
# - add the src directory (where the actual code is) to the build

# ----------------------------------------------------------------------------------------#
# General configuration
# ----------------------------------------------------------------------------------------#

# version 3.18 introduces imported target LAPACK::LAPACK
cmake_minimum_required(VERSION 3.18...3.25)

project(mappraiser
        VERSION 2.0
        DESCRIPTION "MAPPRAISER: A massively parallel map-making framework for multi-kilo pixel CMB experiments"
        LANGUAGES C CXX)

# Enforce language standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# We are building libraries that will eventually be linked into shared modules.
# All code should be built with PIC.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# ----- RPATH -----

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

# ----------------------------------------------------------------------------------------#
# External packages
# ----------------------------------------------------------------------------------------#

find_package(MIDAPACK)

find_package(MPI REQUIRED)

find_package(CFITSIO)

option(ECG "Build mappraiser with ECG solver" OFF)
option(USE_MKL "Use MKL as provider for BLAS/LAPACK libraries" OFF)

if (ECG)
    # ECG solver requires the use of MKL
    set(USE_MKL ON)

    # preAlps + metis dependencies
    # TODO simplify and move to a cmake module

    # create interface library with ECG dependencies
    add_library(ecg_deps INTERFACE)

    if (NOT DEFINED ENV{PREALPSROOT})

        message(FATAL_ERROR "ECG solver needs preAlps. The environment variable PREALPSROOT must be set.")

    elseif (NOT DEFINED ENV{METISROOT})

        message(FATAL_ERROR "ECG solver depends on metis. The environment variable METISROOT must be set.")

    else ()
        # all good, paths to preAlps and metis are provided

        # components of preAlps as imported static libraries
        add_library(preAlps_release STATIC IMPORTED)
        add_library(cplm_core STATIC IMPORTED)
        add_library(cplm_light STATIC IMPORTED)
        add_library(cplm_v0 STATIC IMPORTED)

        set_property(
                TARGET preAlps_release
                PROPERTY IMPORTED_LOCATION $ENV{PREALPSROOT}/lib/libpreAlps_release.a
        )
        set_property(
                TARGET cplm_core
                PROPERTY IMPORTED_LOCATION $ENV{PREALPSROOT}/lib/libcplm_core.a
        )
        set_property(
                TARGET cplm_light
                PROPERTY IMPORTED_LOCATION $ENV{PREALPSROOT}/lib/libcplm_light.a
        )
        set_property(
                TARGET cplm_v0
                PROPERTY IMPORTED_LOCATION $ENV{PREALPSROOT}/lib/libcplm_v0.a
        )

        # metis as imported shared library
        add_library(metis SHARED IMPORTED)
        set_property(TARGET metis PROPERTY IMPORTED_LOCATION $ENV{METISROOT}/lib/libmetis.so)

        # link using the interface library
        target_compile_options(ecg_deps
                INTERFACE
                -I$ENV{PREALPSROOT}/utils/cplm_core
                -I$ENV{PREALPSROOT}/utils/cplm_light
                -I$ENV{PREALPSROOT}/utils/cplm_v0
                -I$ENV{PREALPSROOT}/utils
                -I$ENV{PREALPSROOT}/src/preconditioners
                -I$ENV{PREALPSROOT}/src/solvers
                -I$ENV{METISROOT}/include
                )
        target_link_libraries(ecg_deps INTERFACE preAlps_release cplm_light cplm_v0 cplm_core metis)

    endif ()

endif () # ECG

# First look for MKL, since it will provide support for BLAS/LAPACK
if (USE_MKL)
    find_package(MKL REQUIRED)
endif ()

if (MKL_FOUND)
    # Use MKL for BLAS / LAPACK
    set(BLAS_LIBRARIES "${MKL_LIBRARIES}")
    set(LAPACK_LIBRARIES "${MKL_LIBRARIES}")
    set(BLAS_FOUND TRUE)
    set(LAPACK_FOUND TRUE)
else ()
    # We need a BLAS/LAPACK implementation other than MKL
    find_package(BLAS)
    find_package(LAPACK REQUIRED)
endif ()

# ----------------------------------------------------------------------------------------#
# Subdirectories
# ----------------------------------------------------------------------------------------#

add_subdirectory(src)

# TODO handle Python package in a dedicated CMakeLists file?
# add_subdirectory(python)
