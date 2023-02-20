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

if(USE_ECG)
    set(SRC_LIST APPEND "ecg.c")
endif()

add_library(mappraiser SHARED ${SRC_LIST} ${HEADER_LIST})

find_package(LAPACK REQUIRED)

target_include_directories(
    mappraiser 
    PRIVATE ${CMAKE_SOURCE_DIR}/include/midapack
            ${CMAKE_SOURCE_DIR}/include/mappraiser
)

target_link_libraries(
    mappraiser
    PRIVATE mapp_compile_options
    PUBLIC midapack
           MPI::MPI_C
           cfitsio
#           lapack
)

# All users of this library will need C99
target_compile_features(mappraiser PUBLIC c_std_99)

# IDEs should put the headers in a nice place
source_group(
    TREE "${PROJECT_SOURCE_DIR}/include"
    PREFIX "mappraiser Header Files"
    FILES ${HEADER_LIST}
)
