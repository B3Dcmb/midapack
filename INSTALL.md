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
cmake -S . -B build --install-prefix <path>
cmake --build build
cmake --install build
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

It can be installed from the mappraiser subdirectory, *once midapack is installed*:

```
cd mappraiser
cmake -S . -B build --install-prefix <path>
cmake --build build
cmake --install build
```

### Dependencies

Mappraiser requires the following libraries:

- MIDAPACK
- MPI
- LAPACK
- CFITSIO

For the installation to work properly, make sure that the environment variable `MIDAPACKROOT` is set to whatever installation prefix was used when installing midapack.

The user may want to use a LAPACK implementation provided by Intel MKL (Math Kernel Library).
If so, the feature may be enabled by passing the option `-D USE_MKL=ON`.

### ECG solver

The ECG (enlarged conjugate gradient) solver is optionally built by passing the option `-D ECG=ON` to cmake.

In that case, Mappraiser will also need:

- preAlps (https://github.com/NLAFET/preAlps)
- METIS

Both libraries' locations are to be specified through environment variables `<library>ROOT`, which can be set for a one-time use:

```
PREALPSROOT=<path> METISROOT=<path> cmake [...]
```

In that case, since preAlps is compiled with Intel MKL, Mappraiser will search for an
MKL-provided LAPACK implementation.