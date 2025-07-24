#.rst:
# FindKokkos
# -------
#
# Find Kokkos library
#
# This module finds if Kokkos is installed and selects a default
# configuration for use, preference given to a shared library if both
# shared and static libraries are found
#
# The following variables are optionally searched for defaults
#  KOKKOS_ROOT_DIR:    Base directory where all Kokkos components are found
#
# The following are set after configuration is done:
#  KOKKOS_FOUND
#  KOKKOS_INCLUDE_DIRS
#  KOKKOS_LIBRARIES
#  KOKKOS_LIBRARY_DIRS

find_package(PkgConfig QUIET)

# Check for Kokkos installed via pkg-config
if(PKG_CONFIG_FOUND)
  pkg_check_modules(PC_KOKKOS QUIET kokkos)
endif()

# Set default search paths
set(KOKKOS_ROOT_PATHS
  ${KOKKOS_ROOT_DIR}
  $ENV{KOKKOS_ROOT}
  /usr/local
  /usr
)

# Find the headers
find_path(KOKKOS_INCLUDE_DIR
  NAMES Kokkos_Core.hpp
  HINTS ${PC_KOKKOS_INCLUDE_DIRS}
  PATHS ${KOKKOS_ROOT_PATHS}
  PATH_SUFFIXES include
)

# Find the library
find_library(KOKKOS_LIBRARY
  NAMES kokkoscore
  HINTS ${PC_KOKKOS_LIBRARY_DIRS}
  PATHS ${KOKKOS_ROOT_PATHS}
  PATH_SUFFIXES lib lib64
)

# Handle standard args
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Kokkos
  REQUIRED_VARS KOKKOS_LIBRARY KOKKOS_INCLUDE_DIR
  VERSION_VAR PC_KOKKOS_VERSION
)

if(KOKKOS_FOUND)
  set(KOKKOS_INCLUDE_DIRS ${KOKKOS_INCLUDE_DIR})
  set(KOKKOS_LIBRARIES ${KOKKOS_LIBRARY})
  get_filename_component(KOKKOS_LIBRARY_DIRS ${KOKKOS_LIBRARY} DIRECTORY)
  
  # Create imported target
  if(NOT TARGET Kokkos::kokkos)
    add_library(Kokkos::kokkos UNKNOWN IMPORTED)
    set_target_properties(Kokkos::kokkos PROPERTIES
      IMPORTED_LOCATION "${KOKKOS_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${KOKKOS_INCLUDE_DIR}"
    )
  endif()
endif()

mark_as_advanced(KOKKOS_INCLUDE_DIR KOKKOS_LIBRARY)