# LIBMIDAPACK
# - add target for midapack
# - use an interface library to specify compile options and dependencies
# - add build subdirectories for the midapack components (mapmat, toeplitz, templates)
# - install libmidapack on the system

# reference the object files of the midapack components
# using the generator expression $<TARGET_OBJECTS:>
add_library(midapack SHARED $<TARGET_OBJECTS:mapmat> $<TARGET_OBJECTS:toeplitz>)

# TODO build both a shared and static library? allow to choose?

# ----------------------------------------------------------------------------------------#
# create interface target with compile options and dependencies
# ----------------------------------------------------------------------------------------#

add_library(midapack_base INTERFACE)

include(CheckCCompilerFlag)

# macro that checks if flag if supported for C, if so add to midapack_base
macro(ADD_C_FLAG_IF_AVAIL_MIDAPACK FLAG)
    if (NOT "${FLAG}" STREQUAL "")
        # create a variable for checking the flag if supported, e.g.:
        # -fp-model=precise --> c_fp_model_precise
        string(REGEX REPLACE "^-" "c_" FLAG_NAME "${FLAG}")
        string(REPLACE "-" "_" FLAG_NAME "${FLAG_NAME}")
        string(REPLACE " " "_" FLAG_NAME "${FLAG_NAME}")
        string(REPLACE "=" "_" FLAG_NAME "${FLAG_NAME}")

        check_c_compiler_flag("${FLAG}" ${FLAG_NAME})

        if (${FLAG_NAME})
            target_compile_options(midapack_base INTERFACE $<$<COMPILE_LANGUAGE:C>:${FLAG}>)
        endif ()
    endif ()
endmacro()

# macro that checks if list of flags are supported for C
macro(ADD_FLAGS_IF_AVAIL_MIDAPACK)
    foreach (FLAG ${ARGN})
        add_c_flag_if_avail_midapack("${FLAG}")
    endforeach ()
endmacro()

# standard flags for C and C++
add_flags_if_avail_midapack("-W" "-Wall" "-Wextra" "-Wshadow")

if (MPI_FOUND)
    target_compile_definitions(midapack_base INTERFACE W_MPI)
    target_link_libraries(midapack_base INTERFACE MPI::MPI_C)
endif ()

if (OpenMP_FOUND)
    target_compile_definitions(midapack_base INTERFACE W_OPENMP)
    target_link_libraries(midapack_base INTERFACE OpenMP::OpenMP_C) # FIXME
endif ()

# ----------------------------------------------------------------------------------------#
# add source directories for midapack components
# ----------------------------------------------------------------------------------------#

add_subdirectory(src/mapmat)
add_subdirectory(src/toeplitz)

# TODO: add template operations to midapack library
# add_subdirectory(src/templates)

target_link_libraries(midapack PRIVATE mapmat toeplitz)

# ----------------------------------------------------------------------------------------#
# install midapack library in CMAKE_INSTALL_PREFIX
# ----------------------------------------------------------------------------------------#

include(GNUInstallDirs)

# install the library
install(TARGETS midapack
        DESTINATION ${CMAKE_INSTALL_LIBDIR})

# install the headers
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
