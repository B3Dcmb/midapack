# Kokkos Integration in MIDAPACK

This document describes the Kokkos integration in MIDAPACK, which provides CPU and GPU acceleration capabilities for high-performance CMB data analysis.

## Overview

The Kokkos integration allows MIDAPACK to run on both CPU and GPU architectures with minimal code changes. Kokkos provides a performance-portable programming model that automatically optimizes for the target architecture.

## Features

- **Transparent Integration**: Existing MIDAPACK code works unchanged
- **Automatic GPU/CPU Execution**: Kokkos automatically selects the best execution space
- **Memory Management**: Automatic host/device memory management
- **Performance Portable**: Same code runs optimally on CPU, GPU, and other architectures
- **Backward Compatibility**: Full compatibility with existing OpenMP and MPI code

## Building with Kokkos Support

### Prerequisites

- CMake 3.18 or later
- C++17 compatible compiler
- Kokkos library

### Installing Kokkos

```bash
# Download and build Kokkos
git clone https://github.com/kokkos/kokkos.git
cd kokkos
mkdir build && cd build

# For CPU-only with OpenMP
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local \
         -DKokkos_ENABLE_OPENMP=ON \
         -DKokkos_ENABLE_SERIAL=ON \
         -DCMAKE_POSITION_INDEPENDENT_CODE=ON

# For GPU support (CUDA example)
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local \
         -DKokkos_ENABLE_CUDA=ON \
         -DKokkos_ENABLE_OPENMP=ON \
         -DKokkos_ARCH_VOLTA70=ON \
         -DCMAKE_POSITION_INDEPENDENT_CODE=ON

make -j4 && sudo make install
```

### Building MIDAPACK with Kokkos

```bash
mkdir build && cd build
cmake .. -DENABLE_KOKKOS=ON -DDISABLE_OPENMP=OFF
make -j4
```

## Usage

### Automatic Acceleration

When Kokkos is enabled, existing MIDAPACK functions automatically use GPU acceleration:

```c
#include "midapack.h"

int main(int argc, char** argv) {
    // Initialize Kokkos (required for GPU acceleration)
    kokkos_initialize(argc, argv);
    
    // Your existing MIDAPACK code works unchanged
    double *V_data = malloc(n * m * sizeof(double));
    double *TV_data = malloc(n * m * sizeof(double));
    double **V = &V_data;
    double **TV = &TV_data;
    
    // This automatically uses GPU if available
    stmm_simple_basic(V, n, m, T, lambda, TV);
    
    // Cleanup
    kokkos_finalize();
    return 0;
}
```

### Manual Control

For more control over execution:

```c
#ifdef HAVE_KOKKOS
#include "kokkos_wrapper.h"

// Check if Kokkos is available
if (kokkos_is_initialized()) {
    // Use Kokkos-accelerated version
    kokkos_stmm_simple_basic(&V, n, m, T, lambda);
} else {
    // Fall back to CPU version
    stmm_simple_basic(V, n, m, T, lambda, TV);
}
#endif
```

## Performance Considerations

### Memory Layout

Kokkos uses optimal memory layouts for the target architecture:
- **CPU**: Row-major layout for cache efficiency
- **GPU**: Column-major layout for coalesced memory access

### Execution Policies

The implementation uses Kokkos parallel execution policies:
- `Kokkos::parallel_for` for loop parallelization
- `MDRangePolicy` for multi-dimensional loops
- Automatic work scheduling and load balancing

### Memory Management

Kokkos Views provide automatic memory management:
- Host/device memory allocation
- Automatic data transfers
- Reference counting for memory safety

## Supported Kernels

Currently supported Kokkos-accelerated kernels:

### Toeplitz Algebra
- `stmm_simple_basic`: Basic Toeplitz matrix-vector multiplication
- Automatically used in higher-level routines like `stmm_core`

### Planned Extensions
- `scmm_basic`: FFT-based convolution kernels
- `stbmm`: Block Toeplitz operations
- Matrix operations in `mapmat` module

## Testing

### Unit Tests

Run the Kokkos integration tests:

```bash
cd build
./test_kokkos_integration
```

### Example Applications

Try the example demonstrating Kokkos acceleration:

```bash
./example_kokkos_usage
```

### Performance Benchmarks

Run comprehensive performance comparisons:

```bash
# Build the benchmark
cd build
cmake .. -DENABLE_KOKKOS=ON  # or OFF for CPU-only comparison
make

# Run the benchmark
./benchmark_kokkos_performance
```

The benchmark tests multiple problem sizes and provides detailed timing comparisons between CPU-only and Kokkos-accelerated execution.

## Troubleshooting

### Common Issues

1. **Kokkos not found**
   ```
   Error: Could NOT find Kokkos
   ```
   Solution: Ensure Kokkos is installed and `PKG_CONFIG_PATH` includes Kokkos

2. **OpenMP conflicts**
   ```
   Error: You enabled Kokkos OpenMP support without enabling OpenMP in the compiler!
   ```
   Solution: Build with `-DDISABLE_OPENMP=OFF` or rebuild Kokkos with `-DKokkos_ENABLE_OPENMP=OFF`

3. **CUDA runtime errors**
   ```
   Error: CUDA runtime error
   ```
   Solution: Ensure CUDA drivers are installed and GPU is accessible

### Debug Mode

Enable debug output:

```bash
export KOKKOS_PRINT_CONFIGURATION=1
./your_application
```

## Architecture Support

### CPU Architectures
- x86_64 with AVX, AVX2, AVX-512
- ARM64 with NEON
- Power9/Power10

### GPU Architectures  
- NVIDIA GPUs (Kepler, Maxwell, Pascal, Volta, Turing, Ampere)
- AMD GPUs (via HIP backend)
- Intel GPUs (via SYCL backend, experimental)

## Performance Results

Performance benefits from Kokkos acceleration vary significantly depending on:
- Problem size (larger problems benefit more from GPU acceleration)
- Hardware configuration (CPU cores, GPU model, memory bandwidth)
- Data layout and memory access patterns
- Kokkos backend configuration (Serial, OpenMP, CUDA, etc.)

### Running Performance Benchmarks

To obtain actual performance measurements for your system:

```bash
# Build the benchmark suite
cd build
cmake .. -DENABLE_KOKKOS=ON -DDISABLE_MPI=ON
make benchmark_kokkos_performance

# Run comprehensive benchmarks
./benchmark_kokkos_performance

# Or use the automated script
cd ..
./run_benchmark.sh
```

The benchmark will test multiple problem sizes and provide detailed comparisons between:
- CPU-only execution (with OpenMP if available)
- Kokkos-accelerated execution (GPU/CPU depending on configuration)

**Important**: Actual performance results depend entirely on your hardware setup. The benchmark must be run on your target system to get meaningful performance data.

## Contributing

When adding new computational kernels:

1. Implement the kernel in `kokkos_wrapper.cpp`
2. Add C wrapper function in `kokkos_wrapper.h`
3. Integrate into existing functions with `#ifdef HAVE_KOKKOS`
4. Add unit tests for the new functionality
5. Update documentation

## Future Roadmap

- [ ] Support for multiple GPU execution
- [ ] Integration with MPI for distributed GPU computing
- [ ] Advanced memory optimization strategies
- [ ] Performance profiling and optimization tools
- [ ] Python API with GPU support