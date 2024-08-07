# LIBMIDAPACK

add_library(midapack SHARED)

# TODO build both a shared and static library? allow to choose?

# ----------------------------------------------------------------------------------------#
# sources files
# ----------------------------------------------------------------------------------------#

add_subdirectory(mapmat)
add_subdirectory(toeplitz)
add_subdirectory(templates)

# ----------------------------------------------------------------------------------------#
# header files
# ----------------------------------------------------------------------------------------#

target_include_directories(midapack PRIVATE
        "../include"
        "../include/mapmat"
        "../include/toeplitz"
        "../include/templates"
        )

#target_include_directories(midapack PUBLIC
#                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libmidapack/include>
#                           $<INSTALL_INTERFACE:include/midapack>
#                           )

# ----------------------------------------------------------------------------------------#
# dependencies (searched for in top-level CMakeLists.txt)
# ----------------------------------------------------------------------------------------#

# MPI

if (MPI_C_FOUND)
    target_compile_definitions(midapack PRIVATE W_MPI)
    target_link_libraries(midapack PUBLIC MPI::MPI_C)
endif ()

# OpenMP

if (OpenMP_C_FOUND)
    target_link_libraries(midapack PUBLIC OpenMP::OpenMP_C)
endif ()

# FFTW (for toeplitz module)

if (FFTW_DOUBLE_LIB_FOUND)
    target_link_libraries(midapack PUBLIC FFTW::Double)
endif ()

# MKL / LAPACK (for templates module)

if (MKL_FOUND)
    target_compile_definitions(midapack PRIVATE HAVE_MKL)
    target_include_directories(midapack PUBLIC "${MKL_INCLUDE_DIRS}")
    target_link_libraries(midapack PUBLIC "${MKL_LIBRARIES}")
elseif (LAPACK_FOUND)
    target_link_libraries(midapack PUBLIC LAPACK::LAPACK)
endif ()

# ----------------------------------------------------------------------------------------#
# install midapack library in CMAKE_INSTALL_PREFIX
# ----------------------------------------------------------------------------------------#

include(GNUInstallDirs)

# install the library
install(TARGETS midapack
        DESTINATION lib)

# install the headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/algebra/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/midapack)
