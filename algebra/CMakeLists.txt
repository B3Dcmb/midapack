# ________________________________________
# External dependencies for core modules

option(DISABLE_MPI OFF)
option(DISABLE_OPENMP OFF)

# Search for MPI and OpenMP if we want to use them
if (NOT DISABLE_MPI)
    find_package(MPI REQUIRED)
endif ()

if (NOT DISABLE_OPENMP)
    find_package(OpenMP REQUIRED)
endif ()

# Search for FFTW (necessary for toeplitz module)
find_package(FFTW REQUIRED)

# ________________________________________
# Subdirectories

add_subdirectory(src)
