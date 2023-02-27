# Installation

This repo contains two projects: Midapack and Mappraiser.

They are built using cmake (>= 3.18). Try `cmake --help` in case of doubt.

## Midapack

In a generic use, Midapack requires few libraries:

- MPI
- OpenMP
- FFTW3

Ensure these libraries are available on your system.

To build the library and install it at a given location (prefix), execute the commands:

```
cmake -S . -B build_midapack --install-prefix <path>
cmake --build build_midapack
cmake --install build_midapack
```

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

The library can be installed by specifying the source directory from the root of the repository:

```
cmake -S mappraiser -B build_mappraiser --install-prefix <path>
cmake --build build_mappraiser
cmake --install build_mappraiser
```

### Dependencies

Mappraiser requires the following libraries:

- MPI
- LAPACK
- CFITSIO

### ECG solver

The ECG (enlarged conjugate gradient) solver is optionally built by passing the option `-D ECG=ON` to cmake.

In that case, Mappraiser will also need:

- preAlps (https://github.com/NLAFET/preAlps)
- METIS

Both libraries' locations need to be specified with environment variables `<library>ROOT`, for example by executing

```
PREALPSROOT=<path> METISROOT=<path> cmake -S mappraiser [...]
```

In that case, since preAlps is compiled with the Intel Kernel Math Library, Mappraiser will search for an
MKL-provided LAPACK implementation.