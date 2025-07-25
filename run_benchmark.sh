#!/bin/bash
# Automated benchmark script for Kokkos performance evaluation

echo "=== MIDAPACK Kokkos Performance Evaluation ==="
echo ""

BUILD_DIR="build"
BENCHMARK_LOG="benchmark_results.txt"

# Function to check dependencies
check_dependencies() {
    echo "Checking build dependencies..."
    
    # Check for required libraries
    if ! pkg-config --exists fftw3; then
        echo "WARNING: FFTW3 not found. You may need to install libfftw3-dev"
        echo "  On Ubuntu/Debian: sudo apt install libfftw3-dev"
        echo "  On CentOS/RHEL: sudo yum install fftw3-devel"
    fi
    
    if ! pkg-config --exists lapacke; then
        echo "WARNING: LAPACKE not found. You may need to install liblapacke-dev"
        echo "  On Ubuntu/Debian: sudo apt install liblapacke-dev"
    fi
    
    echo ""
}

# Check dependencies first
check_dependencies

# Check if build directory exists
if [ ! -d "$BUILD_DIR" ]; then
    echo "Build directory not found. Creating and configuring..."
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    
    # Try with minimal dependencies first
    echo "Configuring with KOKKOS=ON, MPI=OFF, OpenMP=ON..."
    cmake .. -DENABLE_KOKKOS=ON -DDISABLE_MPI=ON -DDISABLE_OPENMP=OFF
    
    if [ $? -ne 0 ]; then
        echo ""
        echo "CMake configuration failed. This may be due to missing dependencies."
        echo "Please install the required development libraries and try again."
        echo ""
        echo "Required packages (Ubuntu/Debian):"
        echo "  sudo apt install libfftw3-dev liblapacke-dev libblas-dev"
        echo ""
        echo "For Kokkos support, you may also need:"
        echo "  sudo apt install libkokkos-dev  # if available"
        echo "  or build Kokkos from source: https://github.com/kokkos/kokkos"
        exit 1
    fi
    
    echo "Building..."
    make -j$(nproc)
    if [ $? -ne 0 ]; then
        echo "ERROR: Build failed"
        exit 1
    fi
    cd ..
fi

# Check if benchmark executable exists
if [ ! -f "$BUILD_DIR/benchmark_kokkos_performance" ]; then
    echo "Benchmark executable not found. Building..."
    cd "$BUILD_DIR"
    make benchmark_kokkos_performance
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to build benchmark"
        echo "Make sure all dependencies are installed and try again"
        exit 1
    fi
    cd ..
fi

echo "Running performance benchmark..."
echo "Results will be saved to: $BENCHMARK_LOG"
echo ""

# Run the benchmark and save results
cd "$BUILD_DIR"
./benchmark_kokkos_performance | tee "../$BENCHMARK_LOG"
BENCHMARK_EXIT_CODE=$?
cd ..

if [ $BENCHMARK_EXIT_CODE -eq 0 ]; then
    echo ""
    echo "Benchmark completed successfully. Results saved to $BENCHMARK_LOG"
else
    echo ""
    echo "Benchmark failed or encountered errors."
fi

echo ""
echo "To compare different configurations, you can:"
echo "1. Rebuild with ENABLE_KOKKOS=OFF for CPU-only baseline:"
echo "   cd build && cmake .. -DENABLE_KOKKOS=OFF -DDISABLE_MPI=ON && make"
echo "2. Test different Kokkos backends (CUDA, OpenMP, etc.)"
echo "3. Vary problem sizes by modifying benchmark_kokkos_performance.c"
echo "4. Run on different hardware configurations"