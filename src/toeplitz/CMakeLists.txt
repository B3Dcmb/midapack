set(
    HEADER_LIST
    "${midapack_SOURCE_DIR}/include/midapack/toeplitz/toeplitz.h"
)

set(
    SRC_LIST
    "toeplitz.c"
    "toeplitz_seq.c"
    "toeplitz_block.c"
    "toeplitz_gappy.c"
    "toeplitz_nofft.c"
    "toeplitz_rshp.c"
    "toeplitz_params.c"
    "toeplitz_utils.c"
    "toeplitz_wizard.c"
)

add_library(toeplitz OBJECT ${SRC_LIST} ${HEADER_LIST})

target_include_directories(toeplitz PRIVATE ${CMAKE_SOURCE_DIR}/include/midapack/toeplitz)

# TODO: link to FFTW3 (not necessary at NERSC thanks to Cray compiler wrappers)
target_link_libraries(
    toeplitz
    PRIVATE mida_compile_options
            OpenMP::OpenMP_C
    PUBLIC MPI::MPI_C
)

# All users of this library will need C99
target_compile_features(toeplitz PUBLIC c_std_99)

# IDEs should put the headers in a nice place
source_group(
    TREE "${PROJECT_SOURCE_DIR}/include/midapack/toeplitz"
    PREFIX "toeplitz Header Files"
    FILES ${HEADER_LIST}
)
