# MIDAPACK

# Top-level CMakeLists file
# - contains the call to the project() command and some global settings
# - search for external dependencies: MPI, OpenMP, FFTW3
# - add relevant subdirectories to the build

# ----------------------------------------------------------------------------------------#
# General configuration
# ----------------------------------------------------------------------------------------#

# Version 3.18 introduces LAPACK::LAPACK
cmake_minimum_required(VERSION 3.18...3.25)

# Project name and a few useful settings
project(midapack
        VERSION 2.1
        DESCRIPTION "Parallel software tools for high performance CMB data analysis"
        LANGUAGES C)

# Enforce language standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Auxiliary files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# We are building libraries that will eventually be linked into shared modules.
# All code should be built with PIC.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ----------------------------------------------------------------------------------------#
# External packages
# ----------------------------------------------------------------------------------------#

option(DISABLE_MPI OFF)
option(DISABLE_OPENMP ON)

# Search for MPI and OpenMP if we want to use them
if (NOT DISABLE_MPI)
    find_package(MPI REQUIRED)
endif ()

if (NOT DISABLE_OPENMP)
    find_package(OpenMP REQUIRED)
endif ()

# Search for FFTW (necessary for toeplitz module)
find_package(FFTW REQUIRED)

# ----------------------------------------------------------------------------------------#
# Subdirectories
# ----------------------------------------------------------------------------------------#

add_subdirectory(libmidapack/src)

# exclude mappraiser-related targets from the build by default
# users have to specify they want to build mappraiser when calling cmake
add_subdirectory(mappraiser EXCLUDE_FROM_ALL)
