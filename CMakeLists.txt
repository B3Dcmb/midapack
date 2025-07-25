# MIDAPACK

# ----------------------------------------------------------------------------------------#
# General configuration
# ----------------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.18)

# Project name and a few useful settings
project(midapack
        VERSION 2.1
        DESCRIPTION "Parallel software tools for high performance CMB data analysis"
        LANGUAGES C CXX)

# Enforce language standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
if ("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif ("${isSystemDir}" STREQUAL "-1")

# ----- Build Type -----
# Code taken from https://cliutils.gitlab.io/modern-cmake/chapters/features.html

set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
            STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
            "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

# ----------------------------------------------------------------------------------------#
# Optional dependencies
# ----------------------------------------------------------------------------------------#

# Kokkos (for GPU/CPU acceleration)
option(ENABLE_KOKKOS "Enable Kokkos for CPU/GPU acceleration" OFF)
if (ENABLE_KOKKOS)
    find_package(Kokkos REQUIRED CONFIG)
    message(STATUS "Found Kokkos: ${Kokkos_DIR}")
endif()

# ----------------------------------------------------------------------------------------#
# Subdirectories
# ----------------------------------------------------------------------------------------#

add_subdirectory(algebra)
add_subdirectory(mappraiser)

# ----------------------------------------------------------------------------------------#
# Example and benchmark executables
# ----------------------------------------------------------------------------------------#

# Kokkos example usage
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/example_kokkos_usage.c")
    add_executable(example_kokkos_usage example_kokkos_usage.c)
    target_include_directories(example_kokkos_usage PRIVATE algebra/include)
    target_link_libraries(example_kokkos_usage midapack m lapacke)
    if(ENABLE_KOKKOS)
        target_link_libraries(example_kokkos_usage Kokkos::kokkos)
        target_compile_definitions(example_kokkos_usage PRIVATE HAVE_KOKKOS)
    endif()
endif()

# Kokkos integration test
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test_kokkos_integration.c")
    add_executable(test_kokkos_integration test_kokkos_integration.c)
    target_include_directories(test_kokkos_integration PRIVATE algebra/include)
    target_link_libraries(test_kokkos_integration midapack m)
    if(ENABLE_KOKKOS)
        target_link_libraries(test_kokkos_integration Kokkos::kokkos)
        target_compile_definitions(test_kokkos_integration PRIVATE HAVE_KOKKOS)
    endif()
endif()

# Kokkos performance benchmark
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_kokkos_performance.c")
    add_executable(benchmark_kokkos_performance benchmark_kokkos_performance.c)
    target_include_directories(benchmark_kokkos_performance PRIVATE algebra/include)
    target_link_libraries(benchmark_kokkos_performance midapack m)
    if(ENABLE_KOKKOS AND TARGET Kokkos::kokkos)
        target_link_libraries(benchmark_kokkos_performance Kokkos::kokkos)
        target_compile_definitions(benchmark_kokkos_performance PRIVATE HAVE_KOKKOS)
    endif()
endif()
