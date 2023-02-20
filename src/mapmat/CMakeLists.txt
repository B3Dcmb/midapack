# Note that headers are optional, and do not affect add_library, but they will not
# show up in IDEs unless they are listed in add_library.
set(
    HEADER_LIST
    "${midapack_SOURCE_DIR}/include/midapack/mapmat/alm.h"
    "${midapack_SOURCE_DIR}/include/midapack/mapmat/als.h"
    "${midapack_SOURCE_DIR}/include/midapack/mapmat/bitop.h"
    "${midapack_SOURCE_DIR}/include/midapack/mapmat/butterfly.h"
    "${midapack_SOURCE_DIR}/include/midapack/mapmat/cindex.h"
    "${midapack_SOURCE_DIR}/include/midapack/mapmat/csort.h"
    "${midapack_SOURCE_DIR}/include/midapack/mapmat/mapmat.h"
    "${midapack_SOURCE_DIR}/include/midapack/mapmat/mapmatc.h"
    "${midapack_SOURCE_DIR}/include/midapack/mapmat/ring.h"
)

set(
    SRC_LIST
    "alm.c"
    "als.c"
    "bitop.c"
    "butterfly.c"
    "cindex.c"
    "csort.c"
    "mapmat.c"
    "mapmatc.c"
    "ring.c"
)

# Make an object library that will be used in the final libmidapack
# Dependencies will be propagated to libmidapack when linking the object library
add_library(mapmat OBJECT ${SRC_LIST} ${HEADER_LIST})

target_include_directories(mapmat PRIVATE ${CMAKE_SOURCE_DIR}/include/midapack/mapmat)

target_link_libraries(
    mapmat
    PRIVATE mida_compile_options
    PUBLIC MPI::MPI_C
)

# IDEs should put the headers in a nice place
source_group(
    TREE "${PROJECT_SOURCE_DIR}/include/midapack/mapmat"
    PREFIX "mapmat Header Files"
    FILES ${HEADER_LIST}
)
