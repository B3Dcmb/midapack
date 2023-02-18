# Works with 3.14 and tested through 3.25
cmake_minimum_required(VERSION 3.14...3.25)

# Project name and a few useful settings
project(
    midapack
    VERSION 2.1
    DESCRIPTION "Parallel software tools for high performance CMB data analysis"
    LANGUAGES C
)

add_library(
    midapack SHARED $<TARGET_OBJECTS:mapmat>
                    $<TARGET_OBJECTS:toeplitz>
)

# in order to see the symbols from the static libraries
# set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--export-all-symbols")

# ----------------------------------------------------------------------------------------#
# set language standards
# ----------------------------------------------------------------------------------------#

set(CMAKE_C_STANDARD    99 CACHE STRING "C language standard")

option(CMAKE_C_STANDARD_REQUIRED    "Require the C language standard to set"    ON)

option(CMAKE_C_EXTENSIONS    "Enable/disable extensions, e.g. -std=gnu11 vs. -std=c11"    ON)

set(CMAKE_POSITION_INDEPENDENT_CODE True)

# ----------------------------------------------------------------------------------------#
# create interface target with compiler flags
# ----------------------------------------------------------------------------------------#

add_library(mida_compile_options INTERFACE)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

include(CheckCCompilerFlag)

# ----------------------------------------------------------------------------------------#
# macro that checks if flag if supported for C, if so add to mida_compile_options
# ----------------------------------------------------------------------------------------#
macro(ADD_C_FLAG_IF_AVAIL FLAG)
    if(NOT "${FLAG}" STREQUAL "")
        # create a variable for checking the flag if supported, e.g.:
        # -fp-model=precise --> c_fp_model_precise
        string(REGEX REPLACE "^-" "c_" FLAG_NAME "${FLAG}")
        string(REPLACE "-" "_" FLAG_NAME "${FLAG_NAME}")
        string(REPLACE " " "_" FLAG_NAME "${FLAG_NAME}")
        string(REPLACE "=" "_" FLAG_NAME "${FLAG_NAME}")

        check_c_compiler_flag("${FLAG}" ${FLAG_NAME})

        if(${FLAG_NAME})
            target_compile_options(mida_compile_options INTERFACE
                $<$<COMPILE_LANGUAGE:C>:${FLAG}>)
        endif()
    endif()
endmacro()

# ----------------------------------------------------------------------------------------#
# macro that checks if list of flags are supported for C
# ----------------------------------------------------------------------------------------#
macro(ADD_FLAGS_IF_AVAIL)
    foreach(FLAG ${ARGN})
        add_c_flag_if_avail("${FLAG}")
        # add_cxx_flag_if_avail("${FLAG}")
    endforeach()
endmacro()

# ----------------------------------------------------------------------------------------#
# options
# ----------------------------------------------------------------------------------------#
option(USE_AVX512 "Enable AVX-512 architecture flags" OFF)
option(USE_SANITIZER "Enable leak detection" OFF)

# standard flags for C and C++
add_flags_if_avail("-W" "-Wall" "-Wextra" "-Wshadow")

# # OpenMP SIMD-only (supported by GCC)
# add_flags_if_avail("-fopenmp-simd")

# enable runtime leak detection
if(USE_SANITIZER)
    add_flags_if_avail("-fsanitize=leak")

    # emit warnings that this feature is not available
    if(NOT c_fsanitize_leak)
        message(WARNING "Sanitizer is not available for selected C compiler")
    endif()

    if(NOT cxx_fsanitize_leak)
        message(WARNING "Sanitizer is not available for selected C++ compiler")
    endif()
endif()

# check for AVX-512 flags
if(USE_AVX512)
    if(CMAKE_C_COMPILER_ID MATCHES "Intel")
        add_flags_if_avail("-xMIC-AVX512")
    else()
        # these flags are supported by newer GCC versions
        add_flags_if_avail("-mavx512f" "-mavx512pf" "-mavx512er" "-mavx512cd")
    endif()
endif()

option(USE_MPI "Enable MPI features" ON)
option(USE_OPENMP "Enable OpenMP multithreading" ON)

if(USE_MPI)
    target_compile_definitions(mida_compile_options INTERFACE W_MPI)
endif()

if(USE_OPENMP)
    target_compile_definitions(mida_compile_options INTERFACE W_OPENMP)
endif()

# ----------------------------------------------------------------------------------------#
# add project subdirectories
# ----------------------------------------------------------------------------------------#

add_subdirectory(src/mapmat)
add_subdirectory(src/toeplitz)
# add_subdirectory(src/mappraiser)
# add_subdirectory(src/templates)
