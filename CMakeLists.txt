# MIDAPACK

# Top-level CMakeLists file
# - contains the call to the project() command and some global settings
# - search for external dependencies: MPI, OpenMP, FFTW3
# - add relevant subdirectories to the build

# ----------------------------------------------------------------------------------------#
# General configuration
# ----------------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.14...3.25)

# Project name and a few useful settings
project(midapack
        VERSION 2.1
        DESCRIPTION "Parallel software tools for high performance CMB data analysis"
        LANGUAGES C Fortran)

# Enforce language standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Auxiliary files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# We are building libraries that will eventually be linked into shared modules.
# All code should be built with PIC.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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

# option(DISABLE_MPI OFF)
# option(DISABLE_OPENMP ON)
# # option(DISABLE_HARMONICS OFF)
# option(USE_MKL "Use MKL as provider for BLAS/LAPACK libraries" OFF)

# # Search for MPI and OpenMP if we want to use them
# if (NOT DISABLE_MPI)
#     find_package(MPI REQUIRED)
# endif ()

# if (NOT DISABLE_OPENMP)
#     find_package(OpenMP REQUIRED)
# endif ()

# # Search for FFTW (necessary for toeplitz module)
# find_package(FFTW REQUIRED)

# find_package(CFITSIO REQUIRED)

# find_package(HEALPIX REQUIRED)

# find_package(S2HAT REQUIRED)


# ----------------------------------------------------------------------------------------#
# Subdirectories
# ----------------------------------------------------------------------------------------#

add_subdirectory(libmidapack)
add_subdirectory(mappraiser)


option(BUILD_Sampling "Build the Power-Spectrum-Sampling as well" ON)

if (BUILD_Sampling)
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/Power-Spectrum-Sampling/src/mappraiser_wrapper_WF
    DESTINATION python)
endif ()