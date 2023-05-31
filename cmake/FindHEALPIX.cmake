# - Check for the presence of HEALPIX
#
# The following variables are set when HEALPIX is found:
#  HEALPIX_FOUND      = Set to true, if all components of HEALPIX have been found.
#  HEALPIX_INCLUDES   = Include path for the header files of HEALPIX
#  HEALPIX_LIBRARIES  = Link these to use HEALPIX
#  HEALPIX_LFLAGS     = Linker flags (optional)

if (NOT HEALPIX_FOUND)

    if (NOT HEALPIX_DIR)
        set(HEALPIX_DIR ${CMAKE_INSTALL_PREFIX})
    endif (NOT HEALPIX_DIR)

    ##____________________________________________________________________________
    ## Check for the header files

    find_path(HEALPIX_INCLUDES
              NAMES chealpix.h
              HINTS $ENV{HEALPIX_DIR} ${CMAKE_INSTALL_PREFIX}
              PATH_SUFFIXES include healpix include/healpix
              )

    ##____________________________________________________________________________
    ## Check for the library

    find_library(HEALPIX_LIBRARIES chealpix
                 HINTS $ENV{HEALPIX_DIR} ${CMAKE_INSTALL_PREFIX}
                 PATH_SUFFIXES lib
                 )
    ##____________________________________________________________________________
    ## Set linker flag

    set(HEALPIX_LFLAGS "-DHEALPIXDATA=${HEALPIXROOT}share/healpix/")

    ##____________________________________________________________________________
    ## Actions taken when all components have been found

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(HEALPIX DEFAULT_MSG HEALPIX_LIBRARIES HEALPIX_INCLUDES)

    if (HEALPIX_FOUND)
        if (NOT HEALPIX_FIND_QUIETLY)
            message(STATUS "Found components for Healpix")
            message(STATUS "HEALPIX_DIR  = ${HEALPIX_DIR}")
            message(STATUS "HEALPIX_INCLUDES  = ${HEALPIX_INCLUDES}")
            message(STATUS "HEALPIX_LIBRARIES = ${HEALPIX_LIBRARIES}")
        endif (NOT HEALPIX_FIND_QUIETLY)
    else (HEALPIX_FOUND)
        if (HEALPIX_FIND_REQUIRED)
            message(FATAL_ERROR "Could not find Healpix!")
        endif (HEALPIX_FIND_REQUIRED)
    endif (HEALPIX_FOUND)



    ##____________________________________________________________________________
    ## Mark advanced variables

    mark_as_advanced(
            HEALPIX_DIR
            HEALPIX_INCLUDES
            HEALPIX_LIBRARIES
    )

endif (NOT HEALPIX_FOUND)