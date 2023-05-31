# - Check for the presence of S2HAT
#
# The following variables are set when S2HAT is found:
#  S2HAT_FOUND      = Set to true, if all components of S2HAT have been found.
#  S2HAT_INCLUDES   = Include path for the header files of S2HAT
#  S2HAT_LIBRARIES  = Link these to use S2HAT
#  S2HAT_LFLAGS     = Linker flags (optional)

if (NOT S2HAT_FOUND)

    if (NOT S2HAT_DIR)
        set(S2HAT_DIR ${CMAKE_INSTALL_PREFIX})
    endif (NOT S2HAT_DIR)

    ##____________________________________________________________________________
    ## Check for the header files

    find_path(S2HAT_INCLUDES
              NAMES s2hat.h
              HINTS $ENV{S2HAT_DIR} ${CMAKE_INSTALL_PREFIX}
              PATH_SUFFIXES include s2hat include/s2hat
              )

    ##____________________________________________________________________________
    ## Check for the library

    find_library(S2HAT_LIBRARIES s2hat_std
                 HINTS $ENV{S2HAT_DIR} ${CMAKE_INSTALL_PREFIX}
                 PATH_SUFFIXES lib
                 )
    ##____________________________________________________________________________
    ## Set linker flag

    set(S2HAT_LFLAGS "-ls2hat_std")
    ##____________________________________________________________________________
    ## Actions taken when all components have been found

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(S2HAT DEFAULT_MSG S2HAT_LIBRARIES S2HAT_INCLUDES)

    if (S2HAT_FOUND)
        if (NOT S2HAT_FIND_QUIETLY)
            message(STATUS "Found components for S2HAT")
            message(STATUS "S2HAT_DIR  = ${S2HAT_DIR}")
            message(STATUS "S2HAT_INCLUDES  = ${S2HAT_INCLUDES}")
            message(STATUS "S2HAT_LIBRARIES = ${S2HAT_LIBRARIES}")
        endif (NOT S2HAT_FIND_QUIETLY)
    else (S2HAT_FOUND)
        if (S2HAT_FIND_REQUIRED)
            message(FATAL_ERROR "Could not find S2HAT!")
        endif (S2HAT_FIND_REQUIRED)
    endif (S2HAT_FOUND)



    ##____________________________________________________________________________
    ## Mark advanced variables

    mark_as_advanced(
            S2HAT_DIR
            S2HAT_INCLUDES
            S2HAT_LIBRARIES
    )

endif (NOT S2HAT_FOUND)