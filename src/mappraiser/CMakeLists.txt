set(
    HEADER_LIST
    "${midapack_SOURCE_DIR}/include/mappraiser/mappraiser.h"
    "${midapack_SOURCE_DIR}/include/midapack/midapack.h"
)

set(
    SRC_LIST
    "createToeplitz.c"
    "iofiles.c"
    "mappraiser.c"
    "pcg_true.c"
    "precond.c"
)

if(ECG)
    list(APPEND SRC_LIST "ecg.c")
    if(NOT DEFINED ENV{PREALPSROOT})
        message(FATAL_ERROR "ECG solver needs preAlps. The environment variable PREALPSROOT must be set.")
    elseif(NOT DEFINED ENV{METISROOT})
        message(FATAL_ERROR "ECG solver depends on metis. The environment variable METISROOT must be set.")
    else()
        # all good, paths to preAlps and metis are provided

        # import static libraries of preAlps
        add_library(preAlps_release STATIC IMPORTED)
        add_library(cplm_core STATIC IMPORTED)
        add_library(cplm_light STATIC IMPORTED)
        add_library(cplm_v0 STATIC IMPORTED)
        set_property(TARGET preAlps_release PROPERTY
                     IMPORTED_LOCATION $ENV{PREALPSROOT}/lib/libpreAlps_release.a)
        set_property(TARGET cplm_core PROPERTY
                     IMPORTED_LOCATION $ENV{PREALPSROOT}/lib/libcplm_core.a)
        set_property(TARGET cplm_light PROPERTY
                     IMPORTED_LOCATION $ENV{PREALPSROOT}/lib/libcplm_light.a)
        set_property(TARGET cplm_v0 PROPERTY
                     IMPORTED_LOCATION $ENV{PREALPSROOT}/lib/libcplm_v0.a)

        # same for metis
        add_library(metis SHARED IMPORTED)
        set_property(TARGET metis PROPERTY IMPORTED_LOCATION $ENV{METISROOT}/lib/libmetis.so)

        target_compile_options(
            mappraiser_base
            INTERFACE -I$ENV{PREALPSROOT}/utils/cplm_core
                      -I$ENV{PREALPSROOT}/utils/cplm_light
                      -I$ENV{PREALPSROOT}/utils/cplm_v0
                      -I$ENV{PREALPSROOT}/utils
                      -I$ENV{PREALPSROOT}/src/preconditioners
                      -I$ENV{PREALPSROOT}/src/solvers
                      -I$ENV{METISROOT}/include
        )
        target_link_libraries(
            mappraiser_base
            INTERFACE preAlps_release cplm_light cplm_v0 cplm_core metis
        )
    endif()
endif()

if(USE_MKL)
    if(NOT DEFINED ENV{MKLROOT})
        message(FATAL_ERROR "The environment variable MKLROOT must be set in order to use mkl routines.")
    else()
        # path to MKL library is provided
        target_compile_options(
            mappraiser_base
            INTERFACE -DMKL_ILP64 -I$ENV{MKLROOT}/include
        )
        target_link_options(
            mappraiser_base
            INTERFACE -L$ENV{MKLROOT}/lib/intel64 -lmkl_rt
# -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
        )
    endif()
else()
    # do not use MKL, so we look for another LAPACK implementation
    find_package(LAPACK REQUIRED)
    target_link_libraries(mappraiser_base INTERFACE lapack)
endif()

add_library(mappraiser SHARED ${SRC_LIST} ${HEADER_LIST})

target_include_directories(
    mappraiser 
    PRIVATE ${CMAKE_SOURCE_DIR}/include/midapack
            ${CMAKE_SOURCE_DIR}/include/mappraiser
)

target_link_libraries(
    mappraiser
    PRIVATE mappraiser_base
    PUBLIC midapack
           MPI::MPI_C
           cfitsio
)

# All users of this library will need C99
target_compile_features(mappraiser PUBLIC c_std_99)

# IDEs should put the headers in a nice place
source_group(
    TREE "${PROJECT_SOURCE_DIR}/include"
    PREFIX "mappraiser Header Files"
    FILES ${HEADER_LIST}
)
