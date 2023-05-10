# Installation

This repository contains the code of two libraries:

- Midapack's core algebra routines (`libmidapack`)
- Mappraiser (`libmappraiser`)

Both those libraries are built with cmake.
Try `cmake --help` in case of doubt.

## Dependencies

The following libraries are needed:

- MPI
- OpenMP
- FFTW3
- LAPACK (if not using MKL)
- CFITSIO
- MKL (only for ECG solver)

## Building the project

To build the libraries and install them at a given location (prefix),
execute the commands:

```
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=<prefix>
cmake --build build
cmake --install build
```

(cmake >= 3.21 also has the command-line option `--install-prefix`)

It may be useful to prepend the installation path to your `LD_LIBARY_PATH`
by adding the following line in your `.bashrc` (or equivalent):
```
export LD_LIBRARY_PATH="${PREFIX}/lib:${LD_LIBRARY_PATH}"
```

The user may want to use a LAPACK implementation provided by Intel MKL (Math Kernel Library).
If so, the feature may be enabled by passing the option `-D USE_MKL=ON`.

## ECG solver

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

It may be useful to prepend this path to your `PYTHONPATH`
by adding the following line in your `.bashrc` (or equivalent):
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

## examples (deprecated)

To build the examples binaries for the core library, do:

- make mapmat_example to generate the examples from the mapmat module
- make toeplitz_example to generate the examples from the Toeplitz module

