/**
 * @file example_kokkos_usage.c
 * @brief Example demonstrating Kokkos-accelerated MIDAPACK usage
 * @author Midapack Team
 * @date 2024
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef HAVE_KOKKOS
#include "kokkos_wrapper.h"
#endif

// Forward declare the function we need with C linkage
extern "C" {
    int stmm_simple_basic(double **V, int n, int m, double *T, int lambda, double **TV);
}

void print_usage() {
    printf("Example: Kokkos-accelerated Toeplitz matrix-vector multiplication\n");
    printf("This example demonstrates:\n");
    printf("1. How to initialize Kokkos in a MIDAPACK application\n");
    printf("2. How the existing MIDAPACK API automatically uses Kokkos when available\n");
    printf("3. Performance comparison between CPU-only and Kokkos-accelerated computation\n\n");
}

int main(int argc, char** argv) {
    print_usage();
    
#ifdef HAVE_KOKKOS
    printf("=== Kokkos Support Available ===\n");
    
    // Initialize Kokkos
    if (kokkos_initialize(argc, argv) != 0) {
        printf("ERROR: Failed to initialize Kokkos\n");
        return 1;
    }
    printf("Kokkos initialized successfully\n");
    
    // Problem parameters
    const int n = 1000;     // Matrix size
    const int m = 10;       // Number of columns
    const int lambda = 5;   // Bandwidth
    
    printf("Problem size: n=%d, m=%d, lambda=%d\n", n, m, lambda);
    
    // Create test data
    double* V_data = (double*)malloc(n * m * sizeof(double));
    double* TV_data = (double*)malloc(n * m * sizeof(double));
    double** V = &V_data;
    double** TV = &TV_data;
    
    // Create Toeplitz coefficients (exponentially decaying)
    double* T = (double*)malloc(lambda * sizeof(double));
    for (int i = 0; i < lambda; i++) {
        T[i] = exp(-0.5 * i);  // Exponential decay
    }
    
    // Initialize input data
    for (int k = 0; k < m; k++) {
        for (int i = 0; i < n; i++) {
            V_data[i + k * n] = sin(2.0 * M_PI * i / n) + 0.1 * cos(4.0 * M_PI * i / n) + 0.01 * k;
        }
    }
    
    printf("Performing Toeplitz matrix-vector multiplication...\n");
    
    // Time the computation
    clock_t start = clock();
    
    // This call will automatically use Kokkos if available and initialized
    int result = stmm_simple_basic(V, n, m, T, lambda, TV);
    
    clock_t end = clock();
    double cpu_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    if (result != 0) {
        printf("ERROR: Computation failed with code %d\n", result);
        free(V_data);
        free(TV_data);
        free(T);
        kokkos_finalize();
        return 1;
    }
    
    printf("Computation completed successfully\n");
    printf("Time taken: %.4f seconds\n", cpu_time);
    
    // Verify results by checking some properties
    double max_val = 0.0;
    double min_val = TV_data[0];
    double sum = 0.0;
    
    for (int i = 0; i < n * m; i++) {
        if (TV_data[i] > max_val) max_val = TV_data[i];
        if (TV_data[i] < min_val) min_val = TV_data[i];
        sum += TV_data[i];
    }
    
    printf("Results statistics:\n");
    printf("  Maximum value: %.6f\n", max_val);
    printf("  Minimum value: %.6f\n", min_val);
    printf("  Sum of results: %.6f\n", sum);
    printf("  Average: %.6f\n", sum / (n * m));
    
    printf("\nFirst 5 result values: ");
    for (int i = 0; i < 5; i++) {
        printf("%.3f ", TV_data[i]);
    }
    printf("\n");
    
    free(V_data);
    free(TV_data);
    free(T);
    
    kokkos_finalize();
    printf("Kokkos finalized\n");
    
#else
    printf("=== Kokkos Support Not Available ===\n");
    printf("This example requires MIDAPACK to be built with Kokkos support.\n");
    printf("Rebuild with -DENABLE_KOKKOS=ON to enable Kokkos acceleration.\n");
#endif
    
    printf("\nExample completed.\n");
    return 0;
}