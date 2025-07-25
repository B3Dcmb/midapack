/**
 * @file benchmark_kokkos_performance.c
 * @brief Performance benchmark for Kokkos-accelerated MIDAPACK operations
 * @author Midapack Team
 * @date 2024
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#ifdef HAVE_KOKKOS
#include "kokkos_wrapper.h"
#endif

// Forward declare the function we need with C linkage
extern "C" {
    int stmm_simple_basic(double **V, int n, int m, double *T, int lambda, double **TV);
}

typedef struct {
    int n;
    int m; 
    int lambda;
    double cpu_time;
    double gpu_time;
    double speedup;
    int success;
} benchmark_result_t;

double get_time_diff(clock_t start, clock_t end) {
    return ((double)(end - start)) / CLOCKS_PER_SEC;
}

void create_test_data(double *V_data, double *T, int n, int m, int lambda) {
    // Initialize input data with a realistic signal
    for (int k = 0; k < m; k++) {
        for (int i = 0; i < n; i++) {
            V_data[i + k * n] = sin(2.0 * M_PI * i / n) + 
                               0.1 * cos(4.0 * M_PI * i / n) + 
                               0.01 * k * sin(6.0 * M_PI * i / n);
        }
    }
    
    // Create Toeplitz coefficients (exponentially decaying)
    for (int i = 0; i < lambda; i++) {
        T[i] = exp(-0.5 * i);
    }
}

int benchmark_cpu_only(int n, int m, int lambda, double *cpu_time) {
    printf("  Running CPU-only benchmark...\n");
    
    // Allocate memory
    double* V_data = (double*)malloc(n * m * sizeof(double));
    double* TV_data = (double*)malloc(n * m * sizeof(double));
    double* T = (double*)malloc(lambda * sizeof(double));
    double** V = &V_data;
    double** TV = &TV_data;
    
    if (!V_data || !TV_data || !T) {
        printf("ERROR: Memory allocation failed\n");
        free(V_data); free(TV_data); free(T);
        return -1;
    }
    
    // Create test data
    create_test_data(V_data, T, n, m, lambda);
    
    // Warm up (single run to avoid cold start effects)
    stmm_simple_basic(V, n, m, T, lambda, TV);
    
    // Benchmark multiple runs for better accuracy
    const int num_runs = 5;
    clock_t total_start = clock();
    
    for (int run = 0; run < num_runs; run++) {
        // Reset input data
        create_test_data(V_data, T, n, m, lambda);
        
        int result = stmm_simple_basic(V, n, m, T, lambda, TV);
        if (result != 0) {
            printf("ERROR: CPU computation failed on run %d\n", run);
            free(V_data); free(TV_data); free(T);
            return -1;
        }
    }
    
    clock_t total_end = clock();
    *cpu_time = get_time_diff(total_start, total_end) / num_runs;
    
    free(V_data);
    free(TV_data);
    free(T);
    
    printf("    CPU time (avg of %d runs): %.6f seconds\n", num_runs, *cpu_time);
    return 0;
}

#ifdef HAVE_KOKKOS
int benchmark_kokkos_accelerated(int n, int m, int lambda, double *gpu_time) {
    printf("  Running Kokkos-accelerated benchmark...\n");
    
    // Allocate memory
    double* V_data = (double*)malloc(n * m * sizeof(double));
    double* TV_data = (double*)malloc(n * m * sizeof(double));
    double* T = (double*)malloc(lambda * sizeof(double));
    double** V = &V_data;
    double** TV = &TV_data;
    
    if (!V_data || !TV_data || !T) {
        printf("ERROR: Memory allocation failed\n");
        free(V_data); free(TV_data); free(T);
        return -1;
    }
    
    // Create test data
    create_test_data(V_data, T, n, m, lambda);
    
    // Warm up (includes device memory allocation and transfer costs)
    stmm_simple_basic(V, n, m, T, lambda, TV);
    
    // Benchmark multiple runs
    const int num_runs = 5;
    clock_t total_start = clock();
    
    for (int run = 0; run < num_runs; run++) {
        // Reset input data
        create_test_data(V_data, T, n, m, lambda);
        
        int result = stmm_simple_basic(V, n, m, T, lambda, TV);
        if (result != 0) {
            printf("ERROR: Kokkos computation failed on run %d\n", run);
            free(V_data); free(TV_data); free(T);
            return -1;
        }
    }
    
    clock_t total_end = clock();
    *gpu_time = get_time_diff(total_start, total_end) / num_runs;
    
    free(V_data);
    free(TV_data);
    free(T);
    
    printf("    Kokkos time (avg of %d runs): %.6f seconds\n", num_runs, *gpu_time);
    return 0;
}
#endif

benchmark_result_t run_benchmark(int n, int m, int lambda) {
    benchmark_result_t result = {n, m, lambda, 0.0, 0.0, 0.0, 0};
    
    printf("Benchmarking problem size: n=%d, m=%d, lambda=%d\n", n, m, lambda);
    
    // Run CPU-only benchmark first
    if (benchmark_cpu_only(n, m, lambda, &result.cpu_time) != 0) {
        printf("CPU benchmark failed\n");
        return result;
    }
    
#ifdef HAVE_KOKKOS
    // Run Kokkos-accelerated benchmark
    if (benchmark_kokkos_accelerated(n, m, lambda, &result.gpu_time) != 0) {
        printf("Kokkos benchmark failed\n");
        return result;
    }
    
    // Calculate speedup
    if (result.gpu_time > 0.0) {
        result.speedup = result.cpu_time / result.gpu_time;
    }
    
    printf("  Speedup: %.2fx\n", result.speedup);
#else
    printf("  Kokkos not available - skipping GPU benchmark\n");
    result.gpu_time = -1.0;
    result.speedup = -1.0;
#endif
    
    result.success = 1;
    printf("\n");
    return result;
}

void print_benchmark_summary(benchmark_result_t *results, int num_results) {
    printf("=== BENCHMARK SUMMARY ===\n\n");
    
    printf("%-12s %-8s %-8s %-12s %-12s %-10s\n", 
           "Problem Size", "m", "lambda", "CPU Time (s)", "GPU Time (s)", "Speedup");
    printf("%-12s %-8s %-8s %-12s %-12s %-10s\n",
           "------------", "----", "------", "-----------", "-----------", "-------");
    
    for (int i = 0; i < num_results; i++) {
        if (!results[i].success) continue;
        
        printf("%-12d %-8d %-8d %-12.6f ", 
               results[i].n, results[i].m, results[i].lambda, results[i].cpu_time);
        
#ifdef HAVE_KOKKOS
        if (results[i].gpu_time > 0) {
            printf("%-12.6f %-10.2fx\n", results[i].gpu_time, results[i].speedup);
        } else {
            printf("%-12s %-10s\n", "N/A", "N/A");
        }
#else
        printf("%-12s %-10s\n", "N/A", "N/A");
#endif
    }
    
    printf("\nNotes:\n");
    printf("- CPU times include OpenMP parallelization when available\n");
    printf("- GPU times include memory transfer overhead\n");
    printf("- Results averaged over multiple runs\n");
    printf("- Performance depends heavily on hardware configuration\n");
}

int main(int argc, char** argv) {
    printf("=== MIDAPACK Kokkos Performance Benchmark ===\n\n");
    
#ifdef HAVE_KOKKOS
    printf("Kokkos support: ENABLED\n");
    
    // Initialize Kokkos
    if (kokkos_initialize(argc, argv) != 0) {
        printf("ERROR: Failed to initialize Kokkos\n");
        return 1;
    }
    printf("Kokkos initialized successfully\n\n");
#else
    printf("Kokkos support: DISABLED\n");
    printf("Note: Only CPU benchmarks will be run\n\n");
#endif
    
    // Define benchmark cases (problem sizes to test)
    struct {
        int n, m, lambda;
    } test_cases[] = {
        {100, 5, 3},      // Small problem
        {1000, 10, 5},    // Medium problem  
        {5000, 10, 7},    // Large problem
        {10000, 5, 5},    // Very large problem (if system can handle it)
    };
    
    int num_cases = sizeof(test_cases) / sizeof(test_cases[0]);
    benchmark_result_t *results = (benchmark_result_t*)malloc(num_cases * sizeof(benchmark_result_t));
    
    if (!results) {
        printf("ERROR: Failed to allocate memory for results\n");
#ifdef HAVE_KOKKOS
        kokkos_finalize();
#endif
        return 1;
    }
    
    // Run all benchmarks
    for (int i = 0; i < num_cases; i++) {
        results[i] = run_benchmark(test_cases[i].n, test_cases[i].m, test_cases[i].lambda);
    }
    
    // Print summary
    print_benchmark_summary(results, num_cases);
    
    free(results);
    
#ifdef HAVE_KOKKOS
    kokkos_finalize();
    printf("\nKokkos finalized\n");
#endif
    
    printf("Benchmark completed.\n");
    return 0;
}