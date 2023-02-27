# Installation

This repo contains two projects: Midapack and Mappraiser.

They are built using cmake (version>=3.18 needed).

## Midapack

In a generic use, Midapack requires few libraries: MPI, OpenMP and FFTW3.
Ensure these libraries are available on your system.

To build the library and install it at a given location (prefix), execute the commands:
```
cmake -S . -B build_midapack --install-prefix <path>
cmake --build build_midapack
cmake --install build_midapack
```
The shared library will then be found at `prefix/lib/libmidapack.so` and the headers in `prefix/include/midapack`.

It is possible to disable MPI (respectively OpenMP) by passing the option `-D DISABLE_MPI` (resp. `-D DISABLE_OPENMP`).

[//]: # (TODO Add help text to display when calling cmake --help ?)

### examples (deprecated)

To build the examples binaries for the library, do:
- make mapmat_example  to generate the examples from the mapmat module
- make toeplitz_example to generate the examples from the Toeplitz module

## Mappraiser

Much in the same way, Mappraiser may be installed by executing:
```
cmake -S . -B build_mappraiser --install-prefix <path>
cmake --build build_mappraiser
cmake --install build_mappraiser
```

