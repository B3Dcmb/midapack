# MAPPRAISER

# ----------------------------------------------------------------------------------------#
# External packages
# ----------------------------------------------------------------------------------------#

find_package(MPI REQUIRED)

find_package(CFITSIO REQUIRED)

option(ECG "Build mappraiser with ECG solver" OFF)
option(MKL "Use MKL as provider for BLAS/LAPACK libraries" OFF)

if (ECG)
    # ECG solver requires the use of MKL
    set(MKL ON)

    # preAlps + metis dependencies
    # TODO simplify and move to a cmake module

    find_package(METIS REQUIRED)

    # create interface library with ECG dependencies
    add_library(preAlps INTERFACE)

    if (NOT DEFINED ENV{PREALPS_ROOT})

        message(FATAL_ERROR "ECG solver needs preAlps. The environment variable PREALPS_ROOT must be set.")

    else ()
        # all good, path to preAlps is provided

        # components of preAlps as imported static libraries
        add_library(preAlps_release STATIC IMPORTED)
        add_library(cplm_core STATIC IMPORTED)
        add_library(cplm_light STATIC IMPORTED)
        add_library(cplm_v0 STATIC IMPORTED)

        set_property(
                TARGET preAlps_release
                PROPERTY IMPORTED_LOCATION $ENV{PREALPS_ROOT}/lib/libpreAlps_release.a
        )
        set_property(
                TARGET cplm_core
                PROPERTY IMPORTED_LOCATION $ENV{PREALPS_ROOT}/lib/libcplm_core.a
        )
        set_property(
                TARGET cplm_light
                PROPERTY IMPORTED_LOCATION $ENV{PREALPS_ROOT}/lib/libcplm_light.a
        )
        set_property(
                TARGET cplm_v0
                PROPERTY IMPORTED_LOCATION $ENV{PREALPS_ROOT}/lib/libcplm_v0.a
        )

        # link using the interface library
        target_compile_options(preAlps
                INTERFACE
                -I$ENV{PREALPS_ROOT}/utils/cplm_core
                -I$ENV{PREALPS_ROOT}/utils/cplm_light
                -I$ENV{PREALPS_ROOT}/utils/cplm_v0
                -I$ENV{PREALPS_ROOT}/utils
                -I$ENV{PREALPS_ROOT}/src/preconditioners
                -I$ENV{PREALPS_ROOT}/src/solvers
                )
        target_link_libraries(preAlps INTERFACE preAlps_release cplm_light cplm_v0 cplm_core)

    endif ()

endif () # ECG

# First look for MKL, since it will provide support for BLAS/LAPACK
if (MKL)
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
