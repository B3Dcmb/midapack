# LIBMAPPRAISER
# - add target for mappraiser
# - use an interface library to specify compile options and dependencies
# - add source files and headers to the library
# - install libmappraiser on the system


set(HEADER_LIST
        "${CMAKE_SOURCE_DIR}/libmidapack/include/midapack.h"
        "${PROJECT_SOURCE_DIR}/mappraiser/include/mappraiser.h"
        )

set(SRC_LIST
        "createToeplitz.c"
        "iofiles.c"
        "mappraiser.c"
        "pcg_true.c"
        "precond.c"
        )

add_library(mappraiser SHARED ${SRC_LIST} ${HEADER_LIST})

# ----------------------------------------------------------------------------------------#
# create interface target with compiler flags and dependencies
# ----------------------------------------------------------------------------------------#

add_library(mappraiser_base INTERFACE)

if (MPI_FOUND)
    target_link_libraries(mappraiser_base INTERFACE MPI::MPI_C)
endif ()

if (ECG)
    list(APPEND SRC_LIST "ecg.c")
    target_compile_definitions(mappraiser_base INTERFACE WITH_ECG)
    target_link_libraries(mappraiser_base INTERFACE ecg_deps)
endif ()

if (MKL_FOUND)
    target_compile_definitions(mappraiser_base INTERFACE HAVE_MKL)
endif ()

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

# macros that check if flag if supported for C and C++, if so add to mappraiser_base compile options
macro(ADD_C_FLAG_IF_AVAIL FLAG)
    if (NOT "${FLAG}" STREQUAL "")
        # create a variable for checking the flag if supported, e.g.:
        #   -fp-model=precise --> c_fp_model_precise
        string(REGEX REPLACE "^-" "c_" FLAG_NAME "${FLAG}")
        string(REPLACE "-" "_" FLAG_NAME "${FLAG_NAME}")
        string(REPLACE " " "_" FLAG_NAME "${FLAG_NAME}")
        string(REPLACE "=" "_" FLAG_NAME "${FLAG_NAME}")

        check_c_compiler_flag("${FLAG}" ${FLAG_NAME})

        if (${FLAG_NAME})
            target_compile_options(mappraiser_base INTERFACE $<$<COMPILE_LANGUAGE:C>:${FLAG}>)
        endif ()
    endif ()
endmacro()

macro(ADD_CXX_FLAG_IF_AVAIL FLAG)
    if (NOT "${FLAG}" STREQUAL "")
        # create a variable for checking the flag if supported, e.g.:
        #   -fp-model=precise --> cxx_fp_model_precise
        string(REGEX REPLACE "^-" "cxx_" FLAG_NAME "${FLAG}")
        string(REPLACE "-" "_" FLAG_NAME "${FLAG_NAME}")
        string(REPLACE " " "_" FLAG_NAME "${FLAG_NAME}")
        string(REPLACE "=" "_" FLAG_NAME "${FLAG_NAME}")

        # runs check to see flag is supported by compiler
        check_cxx_compiler_flag("${FLAG}" ${FLAG_NAME})
        if (${FLAG_NAME})
            target_compile_options(mappraiser_base INTERFACE $<$<COMPILE_LANGUAGE:CXX>:${FLAG}>)
        endif ()
    endif ()
endmacro()

# macro that checks if list of flags are supported for C and C++
macro(ADD_FLAGS_IF_AVAIL)
    foreach (FLAG ${ARGN})
        add_c_flag_if_avail("${FLAG}")
        add_cxx_flag_if_avail("${FLAG}")
    endforeach ()
endmacro()

# standard flags for C and C++
add_flags_if_avail("-W" "-Wall" "-Wextra" "-Wshadow")

option(USE_AVX512 "Enable AVX-512 architecture flags" OFF)
option(USE_SANITIZER "Enable leak detection" OFF)

# enable runtime leak detection
if (USE_SANITIZER)
    add_flags_if_avail("-fsanitize=leak")

    # emit warnings that this feature is not available
    if (NOT c_fsanitize_leak)
        message(WARNING "Sanitizer is not available for selected C compiler")
    endif ()

    if (NOT cxx_fsanitize_leak)
        message(WARNING "Sanitizer is not available for selected C++ compiler")
    endif ()
endif ()

# check for AVX-512 flags
if (USE_AVX512)
    if (CMAKE_C_COMPILER_ID MATCHES "Intel")
        add_flags_if_avail("-xMIC-AVX512")
    else ()
        # these flags are supported by newer GCC versions
        add_flags_if_avail("-mavx512f" "-mavx512pf" "-mavx512er" "-mavx512cd")
    endif ()
endif ()

target_include_directories(
        mappraiser PRIVATE
        ${CMAKE_SOURCE_DIR}/libmidapack/include
        ${PROJECT_SOURCE_DIR}/mappraiser/include
)

target_link_libraries(mappraiser mappraiser_base midapack cfitsio)

# IDEs should put the headers in a nice place
#source_group(
#        TREE "${PROJECT_SOURCE_DIR}/include"
#        PREFIX "mappraiser Header Files"
#        FILES ${HEADER_LIST}
#)

# ----------------------------------------------------------------------------------------#
# install the library in the prefix
# ----------------------------------------------------------------------------------------#

include(GNUInstallDirs)

# install the libraries
install(TARGETS mappraiser
        DESTINATION ${CMAKE_INSTALL_LIBDIR})

# install the headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# install the Python files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/python/
        DESTINATION python)
