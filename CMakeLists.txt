# Works with 3.14 and tested through 3.25
cmake_minimum_required(VERSION 3.14...3.25)

# Project name and a few useful settings
project(
    midapack
    VERSION 2.1
    DESCRIPTION "Parallel software tools for high performance CMB data analysis"
    LANGUAGES C
)

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Optionally set things like CMAKE_CXX_STANDARD, CMAKE_POSITION_INDEPENDENT_CODE here
    set(CMAKE_C_STANDARD 11)
    set(CMAKE_C_STANDARD_REQUIRED True)
    set(CMAKE_POSITION_INDEPENDENT_CODE True)

    # # Let's ensure -std=c++xx instead of -std=g++xx
    # set(CMAKE_CXX_EXTENSIONS OFF)

    # Let's nicely support folders in IDEs
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # # Testing only available if this is the main app
    # # Note this needs to be done in the main CMakeLists
    # # since it calls enable_testing, which must be in the
    # # main CMakeLists.
    # include(CTest)

    # # Docs only available if this is the main app
    # find_package(Doxygen)
    # if(Doxygen_FOUND)
    # add_subdirectory(docs)
    # else()
    # message(STATUS "Doxygen not found, not building docs")
    # endif()
endif()

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

# This is where the source code is
add_subdirectory(src/mapmat)
# add_subdirectory(src/toeplitz)
# add_subdirectory(src/mappraiser)
# add_subdirectory(src/templates)

############################################################

# # add_subdirectory( src/toeplitz )
# # add_subdirectory( src/mappraiser )
# list(APPEND EXTRA_LIBS mapmat)

# # list( APPEND EXTRA_LIBS toeplitz )
# # list( APPEND EXTRA_LIBS mappraiser )

# # add midapack lib
# add_library(midapack SHARED)

# # add binary tree to the search path for include files
# target_include_directories(midapack PUBLIC ${PROJECT_BINARY_DIR})
# target_link_libraries(midapack PUBLIC ${EXTRA_LIBS})

############################################################

# # The executable code is here
# add_subdirectory(apps)

# # Testing only available if this is the main app
# # Emergency override MODERN_CMAKE_BUILD_TESTING provided as well
# if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING)
#     AND BUILD_TESTING)
#     add_subdirectory(tests)
# endif()
