# ________________________________________
# External dependencies for core modules

option(DISABLE_MPI "Compile without MPI enabled" OFF)
option(DISABLE_OPENMP "Compile without OpenMP enabled" ON)
# option(DISABLE_HARMONICS OFF)
option(USE_MKL "Use MKL as provider for BLAS/LAPACK libraries" OFF)

# Search for MPI and OpenMP if we want to use them
if (NOT DISABLE_MPI)
    find_package(MPI REQUIRED)
    find_package(MPI COMPONENTS Fortran REQUIRED)
endif ()

if (NOT DISABLE_OPENMP)
    find_package(OpenMP REQUIRED)
endif ()

# Search for FFTW (necessary for toeplitz module)
find_package(FFTW REQUIRED DOUBLE_LIB)

find_package(CFITSIO REQUIRED)

find_package(HEALPIX REQUIRED)

find_package(S2HAT REQUIRED)


# if (NOT DISABLE_HARMONICS)
#     # search s2hat ...
# endif ()

# ----------------------------------------------------------------------------------------#
# Subdirectories
# ----------------------------------------------------------------------------------------#

add_subdirectory(src)
