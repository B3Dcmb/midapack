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
find_package(FFTW REQUIRED DOUBLE_LIB)

# If MKL is enabled look for MKL, and use it to provide support for BLAS/LAPACK
if (MKL)
    find_package(MKL REQUIRED)
endif ()

if (MKL_FOUND)
    # Use MKL for BLAS / LAPACK
    set(BLAS_LIBRARIES "${MKL_LIBRARIES}")
    set(LAPACK_LIBRARIES "${MKL_LIBRARIES}")
    set(BLAS_FOUND TRUE)
    set(LAPACK_FOUND TRUE)
else ()
    # We need a BLAS/LAPACK implementation other than MKL
    find_package(BLAS)

    find_package(LAPACK REQUIRED)
endif ()

# ________________________________________
# Subdirectories

add_subdirectory(src)
