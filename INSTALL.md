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
- MKL (only required by the ECG solver)

## Building the project

To build the libraries and install them at a given location (prefix),
execute the commands:

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=<prefix>
cmake --build build
cmake --install build
```

(cmake >= 3.21 also has the command-line option `--install-prefix`)

Usually it is useful to specify the MPI C compiler using the option `-DCMAKE_C_COMPILER=...` when calling CMake.

It may be useful to prepend the installation path to your `LD_LIBARY_PATH` by adding the following line in your `.bashrc` (or equivalent):

```bash
export LD_LIBRARY_PATH="<prefix>/lib:${LD_LIBRARY_PATH}"
```

The user may want to use a LAPACK implementation provided by Intel MKL (Math Kernel Library).
If so, the feature may be enabled by passing the option `-D MKL=ON`.

## ECG solver

The ECG (enlarged conjugate gradient) solver is optionally built by passing the option `-D ECG=ON` to cmake.

In that case, Mappraiser will also need:

- MKL (instead of any other LAPACK implementation)
- METIS
- preAlps (<https://github.com/NLAFET/preAlps>)

The location of the preAlps libraries are to be specified through a variable `PREALPS_ROOT`, which may be an environment variable or simply set for a one-time use:

```bash
PREALPS_ROOT=<path> cmake [...]
```

### Python

The mappraiser library comes with a Python wrapper and an interface to the TOAST library (<https://github.com/hpc4cmb/toast>).
These files are installed in `<prefix>/lib/python3.xx/site-packages/pymappraiser` (adapt with the version of Python you are using).

It may be useful to prepend this path to your `PYTHONPATH` by adding the following line in your `.bashrc` (or equivalent):

```bash
export PYTHONPATH="<prefix>/lib/python3.xx/site-packages:${PYTHONPATH}"
```

## Testing the installation

If TOAST is installed on your system, you may test your mappraiser installation by executing the following script:

```python
import toast.mpi
from pymappraiser.mappraiser_test import MappraiserTest

def main():
    mt = MappraiserTest()
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
