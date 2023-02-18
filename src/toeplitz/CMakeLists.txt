set(
    HEADER_LIST
    "${midapack_SOURCE_DIR}/include/toeplitz/toeplitz.h"
)

add_library(
    toeplitz
    OBJECT toeplitz.c
           toeplitz_seq.c
           toeplitz_block.c
           toeplitz_gappy.c
           toeplitz_nofft.c
           toeplitz_rshp.c
           toeplitz_params.c
           toeplitz_utils.c
           toeplitz_wizard.c
    ${HEADER_LIST}
)

# We need this directory, and users of our library will need it too
target_include_directories(toeplitz PUBLIC ../../include/toeplitz)

# This depends on (header only) MPI and OpenMP
target_link_libraries(toeplitz PRIVATE MPI::MPI_C OpenMP::OpenMP_C mida_compile_options)

# All users of this library will need C99
target_compile_features(toeplitz PUBLIC c_std_99)

# IDEs should put the headers in a nice place
source_group(
    TREE "${PROJECT_SOURCE_DIR}/include/toeplitz"
    PREFIX "toeplitz Header Files"
    FILES ${HEADER_LIST}
)
