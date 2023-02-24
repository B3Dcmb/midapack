# TOEPLITZ
# - object library (only compile the source files)
# - will be linked by the midapack target

set(toeplitz_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include/toeplitz")

set(HEADER_LIST "${toeplitz_INCLUDE_DIR}/toeplitz.h")

set(SRC_LIST
        "toeplitz.c"
        "toeplitz_seq.c"
        "toeplitz_block.c"
        "toeplitz_gappy.c"
        "toeplitz_nofft.c"
        "toeplitz_utils.c"
        "toeplitz_wizard.c")

add_library(toeplitz OBJECT ${SRC_LIST} ${HEADER_LIST})

target_include_directories(toeplitz PRIVATE ${toeplitz_INCLUDE_DIR})

# TODO link to FFTW3 (not necessary at NERSC thanks to Cray compiler wrappers)
target_link_libraries(toeplitz midapack_base)

# IDEs should put the headers in a nice place
#source_group(
#        TREE "../../include/toeplitz"
#        PREFIX "toeplitz Header Files"
#        FILES ${HEADER_LIST}
#)
