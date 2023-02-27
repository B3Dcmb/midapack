# - Check for the presence of MIDAPACK
#
# The following variables are set when MIDAPACK is found:
#  MIDAPACK_FOUND      = Set to true, if all components of MIDAPACK have been found.
#  MIDAPACK_INCLUDES   = Include path for the header files of MIDAPACK
#  MIDAPACK_LIBRARIES  = Link these to use MIDAPACK
#  MIDAPACK_LFLAGS     = Linker flags (optional)

if (NOT MIDAPACK_FOUND)

    if (NOT MIDAPACK_DIR)
        set(MIDAPACK_DIR ${CMAKE_INSTALL_PREFIX})
    endif (NOT MIDAPACK_DIR)

    ##____________________________________________________________________________
    ## Check for the header files

    find_path(MIDAPACK_INCLUDES
              NAMES midapack.h
              HINTS ${MIDAPACK_DIR} ${CMAKE_INSTALL_PREFIX}
              PATH_SUFFIXES include midapack include/midapack
              )

    ##____________________________________________________________________________
    ## Check for the library

    find_library(MIDAPACK_LIBRARIES midapack
                 HINTS ${MIDAPACK_DIR} ${CMAKE_INSTALL_PREFIX}
                 PATH_SUFFIXES lib
                 )

    ##____________________________________________________________________________
    ## Actions taken when all components have been found

    include(FindPackageHandleStandardArgs)

    find_package_handle_standard_args(MIDAPACK DEFAULT_MSG MIDAPACK_LIBRARIES MIDAPACK_INCLUDES)

    if (MIDAPACK_FOUND)
        if (NOT MIDAPACK_FIND_QUIETLY)
            message(STATUS "Found components for Midapack")
            message(STATUS "MIDAPACK_DIR  = ${MIDAPACK_DIR}")
            message(STATUS "MIDAPACK_INCLUDES  = ${MIDAPACK_INCLUDES}")
            message(STATUS "MIDAPACK_LIBRARIES = ${MIDAPACK_LIBRARIES}")
        endif (NOT MIDAPACK_FIND_QUIETLY)
    else (MIDAPACK_FOUND)
        if (MIDAPACK_FIND_REQUIRED)
            message(FATAL_ERROR "Could not find Midapack!")
        endif (MIDAPACK_FIND_REQUIRED)
    endif (MIDAPACK_FOUND)

    ##____________________________________________________________________________
    ## Mark advanced variables

    mark_as_advanced(
            MIDAPACK_DIR
            MIDAPACK_INCLUDES
            MIDAPACK_LIBRARIES
    )

endif (NOT MIDAPACK_FOUND)