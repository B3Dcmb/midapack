# LIBMAPPRAISER

add_library(mappraiser SHARED)

# ----------------------------------------------------------------------------------------#
# sources files
# ----------------------------------------------------------------------------------------#

target_sources(mappraiser PRIVATE
        "create_toeplitz.c"
        "iofiles.c"
        "mappraiser.c"
        "pcg_true.c"
        "precond.c"
        )

if (ECG)
    target_sources(mappraiser PRIVATE "ecg.c")
endif ()

# ----------------------------------------------------------------------------------------#
# header files
# ----------------------------------------------------------------------------------------#

target_include_directories(mappraiser PRIVATE
        "../include"
        "../include/mappraiser"
        )

#target_include_directories(mappraiser PUBLIC
#                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
#                           $<INSTALL_INTERFACE:include/mappraiser>  # <prefix>/include/mappraiser
#                           )

# ----------------------------------------------------------------------------------------#
# dependencies
# ----------------------------------------------------------------------------------------#

# MIDAPACK

target_include_directories(mappraiser PRIVATE
        "${CMAKE_SOURCE_DIR}/algebra/include"
        )
target_link_libraries(mappraiser PUBLIC midapack)

# MPI

if (MPI_FOUND)
    target_compile_definitions(mappraiser PRIVATE W_MPI)
    target_link_libraries(mappraiser PUBLIC MPI::MPI_C)
endif ()

# CFITSIO

if (CFITSIO_FOUND)
    target_include_directories(mappraiser PRIVATE "${CFITSIO_INCLUDES}")
    target_link_libraries(mappraiser PUBLIC "${CFITSIO_LIBRARIES}")
else ()
    message(FATAL_ERROR "Could not find cfitsio. Consider providing a hint (-D CFITSIO_DIR <path>).")
endif ()

# preAlps and metis

if (ECG)
    target_compile_definitions(mappraiser PRIVATE WITH_ECG)
    target_link_libraries(mappraiser PRIVATE ecg_deps)
endif ()

# MKL

if (MKL_FOUND)
    target_compile_definitions(mappraiser PRIVATE HAVE_MKL)
    target_include_directories(mappraiser PUBLIC "${MKL_INCLUDE_DIRS}")
    target_link_libraries(mappraiser PUBLIC "${MKL_LIBRARIES}")
elseif (LAPACK_FOUND)
    # BLAS/LAPACK
    #    set_target_properties(mappraiser PROPERTIES LINK_FLAGS
    #                          "${LAPACK_LINKER_FLAGS}" #${BLAS_LINKER_FLAGS}"
    #                          )
    #    target_link_libraries(mappraiser PUBLIC "${LAPACK_LIBRARIES}")
    #    target_link_libraries(mappraiser "${BLAS_LIBRARIES}")
    target_link_libraries(mappraiser PUBLIC LAPACK::LAPACK)
endif ()

# ----------------------------------------------------------------------------------------#
# compile options
# ----------------------------------------------------------------------------------------#

include(CheckCCompilerFlag)

# macros that check if flag if supported for C, if so add to mappraiser compile options
macro(add_c_flag_if_avail FLAG)
    if (NOT "${FLAG}" STREQUAL "")
        # create a variable for checking the flag if supported, e.g.:
        #   -fp-model=precise --> c_fp_model_precise
        string(REGEX REPLACE "^-" "c_" FLAG_NAME "${FLAG}")
        string(REPLACE "-" "_" FLAG_NAME "${FLAG_NAME}")
        string(REPLACE " " "_" FLAG_NAME "${FLAG_NAME}")
        string(REPLACE "=" "_" FLAG_NAME "${FLAG_NAME}")

        check_c_compiler_flag("${FLAG}" ${FLAG_NAME})

        if (${FLAG_NAME})
            target_compile_options(mappraiser PRIVATE $<$<COMPILE_LANGUAGE:C>:${FLAG}>)
        endif ()
    endif ()
endmacro()

# macro that checks if list of flags are supported for C
macro(add_flags_if_avail)
    foreach (FLAG ${ARGN})
        add_c_flag_if_avail("${FLAG}")
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

# ----------------------------------------------------------------------------------------#
# install the library in the prefix
# ----------------------------------------------------------------------------------------#

include(GNUInstallDirs)

# install the libraries
install(TARGETS mappraiser
        DESTINATION lib)
#DESTINATION ${CMAKE_INSTALL_LIBDIR})

# install the headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/mappraiser/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/midapack)

# install the Python files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/mappraiser/python/
        DESTINATION python)
