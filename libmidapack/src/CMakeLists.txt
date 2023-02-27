# LIBMIDAPACK

add_library(midapack SHARED)

# TODO build both a shared and static library? allow to choose?

# ----------------------------------------------------------------------------------------#
# sources files
# ----------------------------------------------------------------------------------------#

add_subdirectory(mapmat)
add_subdirectory(toeplitz)
# add_subdirectory(templates)

# ----------------------------------------------------------------------------------------#
# header files
# ----------------------------------------------------------------------------------------#

target_include_directories(midapack PRIVATE
                           "../include"
                           "../include/mapmat"
                           "../include/toeplitz"
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

if (OpenMP_FOUND)
    target_compile_definitions(midapack PRIVATE W_OPENMP)
    target_link_libraries(midapack PUBLIC OpenMP::OpenMP_C)
endif ()

# FFTW (for toeplitz module)

if (FFTW_FOUND)
    #target_compile_definitions(midapack PRIVATE HAVE_FFTW=1)
    target_include_directories(midapack PRIVATE "${FFTW_INCLUDE_DIRS}")
    target_link_libraries(midapack PUBLIC "${FFTW_LIBRARIES}")
    #    if(FFTW_DOUBLE_THREADS_LIB_FOUND)
    #        target_compile_definitions(midapack PRIVATE HAVE_FFTW_THREADS=1)
    #    endif(FFTW_DOUBLE_THREADS_LIB_FOUND)
endif (FFTW_FOUND)

# ----------------------------------------------------------------------------------------#
# compile options
# ----------------------------------------------------------------------------------------#

include(CheckCCompilerFlag)

# macro that checks if flag if supported for C, if so add to midapack
macro(add_c_flag_if_avail_midapack FLAG)
    if (NOT "${FLAG}" STREQUAL "")
        # create a variable for checking the flag if supported, e.g.:
        # -fp-model=precise --> c_fp_model_precise
        string(REGEX REPLACE "^-" "c_" FLAG_NAME "${FLAG}")
        string(REPLACE "-" "_" FLAG_NAME "${FLAG_NAME}")
        string(REPLACE " " "_" FLAG_NAME "${FLAG_NAME}")
        string(REPLACE "=" "_" FLAG_NAME "${FLAG_NAME}")

        check_c_compiler_flag("${FLAG}" ${FLAG_NAME})

        if (${FLAG_NAME})
            target_compile_options(midapack PRIVATE $<$<COMPILE_LANGUAGE:C>:${FLAG}>)
        endif ()
    endif ()
endmacro()

# macro that checks if list of flags are supported for C
macro(add_flags_if_avail_midapack)
    foreach (FLAG ${ARGN})
        add_c_flag_if_avail_midapack("${FLAG}")
    endforeach ()
endmacro()

# standard flags for C and C++
add_flags_if_avail_midapack("-W" "-Wall" "-Wextra" "-Wshadow")

# ----------------------------------------------------------------------------------------#
# install midapack library in CMAKE_INSTALL_PREFIX
# ----------------------------------------------------------------------------------------#

include(GNUInstallDirs)

# install the library
install(TARGETS midapack
        DESTINATION ${CMAKE_INSTALL_LIBDIR})

# install the headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/libmidapack/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/midapack)
