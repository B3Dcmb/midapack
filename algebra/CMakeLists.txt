# ________________________________________
# External dependencies for core modules

option(DISABLE_MPI "Compile without MPI enabled" OFF)
option(DISABLE_OPENMP "Compile without OpenMP enabled" ON)

# Search for MPI and OpenMP if we want to use them
if (NOT DISABLE_MPI)
    find_package(MPI REQUIRED)
endif ()

if (NOT DISABLE_OPENMP)
    find_package(OpenMP)
endif ()

# Search for FFTW (necessary for toeplitz module)
find_package(FFTW REQUIRED)

# ________________________________________
# Subdirectories

add_subdirectory(src)
