#!/bin/bash
# Simple test script for Kokkos integration

echo "Building test for Kokkos integration..."

# Compile the test
g++ -std=c++17 -O2 -DHAVE_KOKKOS \
    -I/home/runner/work/midapack/midapack/algebra/include \
    -I/usr/local/include \
    -L/home/runner/work/midapack/midapack/build/algebra/src \
    -L/usr/local/lib \
    test_kokkos_integration.c \
    -lmidapack -lkokkoscore -llapacke -llapack -lblas -lfftw3 -ldl -lpthread \
    -o test_kokkos_integration

if [ $? -eq 0 ]; then
    echo "Build successful, running test..."
    LD_LIBRARY_PATH=/home/runner/work/midapack/midapack/build/algebra/src:$LD_LIBRARY_PATH ./test_kokkos_integration
else
    echo "Build failed"
    exit 1
fi