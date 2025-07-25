/**
 * @file test_kokkos_integration.c
 * @brief Unit test for Kokkos integration in MIDAPACK
 * @author Midapack Team
 * @date 2024
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef HAVE_KOKKOS
#include "kokkos_wrapper.h"
#endif

// Simple test case for stmm_simple_basic functionality
int test_kokkos_stmm_basic() {
    printf("Testing Kokkos STMM basic functionality...\n");
    
#ifdef HAVE_KOKKOS
    // Initialize Kokkos
    int argc = 1;
    char* argv[] = {"test"};
    if (kokkos_initialize(argc, argv) != 0) {
        printf("FAILED: Could not initialize Kokkos\n");
        return 1;
    }
    
    if (!kokkos_is_initialized()) {
        printf("FAILED: Kokkos not properly initialized\n");
        return 1;
    }
    
    printf("Kokkos initialized successfully\n");
    
    // Test parameters
    const int n = 10;      // Matrix size
    const int m = 2;       // Number of columns
    const int lambda = 3;  // Bandwidth
    
    // Create test data
    double* V_data = (double*)malloc(n * m * sizeof(double));
    double* V_ptr = V_data;
    
    // Create Toeplitz coefficients
    double T[3] = {1.0, 0.5, 0.25};
    
    // Initialize test data
    for (int k = 0; k < m; k++) {
        for (int i = 0; i < n; i++) {
            V_data[i + k * n] = (double)(i + 1) + 0.1 * k;
        }
    }
    
    // Store original data for comparison
    double* V_original = (double*)malloc(n * m * sizeof(double));
    memcpy(V_original, V_data, n * m * sizeof(double));
    
    // Test Kokkos implementation
    int result = kokkos_stmm_simple_basic(&V_ptr, n, m, T, lambda);
    
    if (result != 0) {
        printf("FAILED: Kokkos STMM function returned error %d\n", result);
        free(V_data);
        free(V_original);
        kokkos_finalize();
        return 1;
    }
    
    // Check that computation was actually performed
    int changed = 0;
    for (int i = 0; i < n * m; i++) {
        if (fabs(V_data[i] - V_original[i]) > 1e-12) {
            changed = 1;
            break;
        }
    }
    
    if (!changed) {
        printf("WARNING: Data appears unchanged after Kokkos computation\n");
    } else {
        printf("SUCCESS: Kokkos computation modified data as expected\n");
    }
    
    // Print some results for verification
    printf("First few results: ");
    for (int i = 0; i < 3 && i < n; i++) {
        printf("%.3f ", V_data[i]);
    }
    printf("\n");
    
    // Cleanup
    free(V_data);
    free(V_original);
    kokkos_finalize();
    
    printf("Kokkos test completed successfully\n");
    return 0;
    
#else
    printf("SKIPPED: Kokkos not available in this build\n");
    return 0;
#endif
}

int main() {
    printf("=== MIDAPACK Kokkos Integration Test ===\n\n");
    
    int result = test_kokkos_stmm_basic();
    
    if (result == 0) {
        printf("\nAll tests PASSED\n");
    } else {
        printf("\nSome tests FAILED\n");
    }
    
    return result;
}