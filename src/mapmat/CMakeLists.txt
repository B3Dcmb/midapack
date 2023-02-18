# Note that headers are optional, and do not affect add_library, but they will not
# show up in IDEs unless they are listed in add_library.

# Optionally glob, but only for CMake 3.12 or later:
# file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${ModernCMakeExample_SOURCE_DIR}/include/mapmat/*.h")
set(
    HEADER_LIST
    "${midapack_SOURCE_DIR}/include/mapmat/alm.h"
    "${midapack_SOURCE_DIR}/include/mapmat/als.h"
    "${midapack_SOURCE_DIR}/include/mapmat/bitop.h"
    "${midapack_SOURCE_DIR}/include/mapmat/butterfly.h"
    "${midapack_SOURCE_DIR}/include/mapmat/cindex.h"
    "${midapack_SOURCE_DIR}/include/mapmat/csort.h"
    "${midapack_SOURCE_DIR}/include/mapmat/mapmat.h"
    "${midapack_SOURCE_DIR}/include/mapmat/mapmatc.h"
    "${midapack_SOURCE_DIR}/include/mapmat/ring.h"
)

# Make an automatic library - will be static or dynamic based on user setting
add_library(mapmat alm.c als.c bitop.c butterfly.c cindex.c csort.c mapmat.c mapmatc.c ring.c ${HEADER_LIST})

# We need this directory, and users of our library will need it too
target_include_directories(mapmat PUBLIC ../../include/mapmat)

# This depends on (header only) MPI and OpenMP
target_link_libraries(mapmat PRIVATE MPI::MPI_C OpenMP::OpenMP_C)

# # All users of this library will need at least C++11
# target_compile_features(mapmat PUBLIC cxx_std_11)

# IDEs should put the headers in a nice place
source_group(
    TREE "${PROJECT_SOURCE_DIR}/include"
    PREFIX "Header Files"
    FILES ${HEADER_LIST}
)
