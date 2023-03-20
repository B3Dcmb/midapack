# Installation

This repo contains two projects: Midapack and Mappraiser.

They are built using cmake. Try `cmake --help` in case of doubt.

## Midapack

In a generic use, Midapack requires few libraries:

- MPI
- OpenMP
- FFTW3

Ensure these libraries are available on your system.

To build the library and install it at a given location (prefix), execute the commands:

```
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=<prefix>
cmake --build build
cmake --install build
```
(cmake >= 3.21 also has the command-line option `--install-prefix`)

The shared library will then be found at `prefix/lib/libmidapack.so` and the headers in `prefix/include/midapack`.

It is possible to disable MPI (respectively OpenMP) by passing the option `-D DISABLE_MPI=ON` (
resp. `-D DISABLE_OPENMP=ON`).

[//]: # (TODO Add help text to display when calling cmake --help ?)

### examples (deprecated)

To build the examples binaries for the library, do:

- make mapmat_example to generate the examples from the mapmat module
- make toeplitz_example to generate the examples from the Toeplitz module

## Mappraiser

Mappraiser has an independent build system, as it will eventually be moved to an other git repository.

It can be installed from the mappraiser subdirectory, *once midapack is installed*:

```
cd mappraiser
cmake -S . -B build --install-prefix=<prefix>
cmake --build build
cmake --install build
```

### Dependencies

Mappraiser requires the following libraries:

- MIDAPACK
- MPI
- LAPACK
- CFITSIO

Typically the path to the Midapack installation can be provided by adding the hint `-D MIDAPACK_DIR=<path>`
to the first cmake command, or by adding `<prefix>/lib` to `$LD_LIBRARY_PATH`.

The user may want to use a LAPACK implementation provided by Intel MKL (Math Kernel Library).
If so, the feature may be enabled by passing the option `-D USE_MKL=ON`.

### ECG solver

The ECG (enlarged conjugate gradient) solver is optionally built by passing the option `-D ECG=ON` to cmake.

In that case, Mappraiser will also need:

- preAlps (https://github.com/NLAFET/preAlps)
- METIS

Both libraries' locations are to be specified through environment variables `<library>ROOT`, which can be set for a
one-time use:

```
PREALPSROOT=<path> METISROOT=<path> cmake [...]
```

In that case, since preAlps is compiled with Intel MKL, Mappraiser will search for an
MKL-provided LAPACK implementation.

### Python

The mappraiser library comes with a Python wrapper and an interface to the TOAST library
(https://github.com/hpc4cmb/toast). These files are installed in `CMAKE_INSTALL_PREFIX/python`.

It may then be useful to prepend this path to your `PYTHONPATH` by having this line in your `.bashrc`
```
export PYTHONPATH="${PREFIX}/python:${PYTHONPATH}"
```

## Testing the installation

If TOAST is installed on your system, you may test your mappraiser installation by executing the following script:
```
import toast.mpi
from TOAST_interface import mappraiser_test

def main():
    mt = mappraiser_test.MappraiserTest()
    mt.setUp()
    mt.test_mappraiser_interface()

if __name__ == "__main__":
    world, procs, rank = toast.mpi.get_world()
    with toast.mpi.exception_guard(comm=world):
        main()
```
