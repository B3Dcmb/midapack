# MAPMAT
# - object library (only compile the source files)
# - will be linked by the midapack target

# Note that headers are optional, and do not affect add_library, but they will not
# show up in IDEs unless they are listed in add_library.
set(mapmat_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include/mapmat")

set(
    HEADER_LIST
    "${mapmat_INCLUDE_DIR}/alm.h"
    "${mapmat_INCLUDE_DIR}/als.h"
    "${mapmat_INCLUDE_DIR}/bitop.h"
    "${mapmat_INCLUDE_DIR}/butterfly.h"
    "${mapmat_INCLUDE_DIR}/cindex.h"
    "${mapmat_INCLUDE_DIR}/csort.h"
    "${mapmat_INCLUDE_DIR}/mapmat.h"
    "${mapmat_INCLUDE_DIR}/mapmatc.h"
    "${mapmat_INCLUDE_DIR}/ring.h"
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

target_include_directories(mapmat PRIVATE ${mapmat_INCLUDE_DIR})

target_link_libraries(mapmat midapack_base)

# IDEs should put the headers in a nice place
#source_group(
#    TREE "../../include/mapmat"
#    PREFIX "mapmat Header Files"
#    FILES ${HEADER_LIST}
#)
