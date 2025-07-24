/**
 * @file kokkos_wrapper.cpp
 * @brief C++ implementation of Kokkos wrapper for C interface
 * @author Midapack Team
 * @date 2024
 */

#include "kokkos_wrapper.h"
#include <Kokkos_Core.hpp>
#include <cstdlib>
#include <cstring>
#include <algorithm>

extern "C" {

int kokkos_initialize(int argc, char** argv) {
    try {
        if (!Kokkos::is_initialized()) {
            Kokkos::initialize(argc, argv);
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

void kokkos_finalize(void) {
    if (Kokkos::is_initialized()) {
        Kokkos::finalize();
    }
}

int kokkos_is_initialized(void) {
    return Kokkos::is_initialized() ? 1 : 0;
}

void* kokkos_malloc(size_t size) {
    try {
        return Kokkos::kokkos_malloc(size);
    } catch (...) {
        return nullptr;
    }
}

void kokkos_free(void* ptr) {
    if (ptr != nullptr) {
        Kokkos::kokkos_free(ptr);
    }
}

void kokkos_memcpy_to_device(void* dst, const void* src, size_t size) {
    // For simplicity, just use memcpy for now
    memcpy(dst, src, size);
}

void kokkos_memcpy_to_host(void* dst, const void* src, size_t size) {
    // For simplicity, just use memcpy for now
    memcpy(dst, src, size);
}

/**
 * @brief Kokkos implementation of Toeplitz matrix-vector multiplication
 * 
 * This function implements the same computation as stmm_simple_basic but using Kokkos
 * for parallel execution on CPU or GPU. The algorithm computes:
 * 
 * For each column k (0 to m-1):
 *   For each row i (offset_edges to n-offset_edges):
 *     TV[i + k*n] = sum(T[|j-i|] * V[j + k*n]) for j in [j_first, j_last)
 * 
 * where j_first = max(i - (lambda-1), 0) and j_last = min(i + lambda, n)
 */
int kokkos_stmm_simple_basic(double** V, int n, int m, const double* T, int lambda) {
    if (!Kokkos::is_initialized()) {
        return -1; // Kokkos not initialized
    }
    
    try {
        // Parameters matching the original algorithm
        const int distcorrmin = lambda - 1;
        const int offset_edges = distcorrmin; // flag_nocomputeedges = 1
        
        // Create Kokkos views using standard allocation
        auto V_view = Kokkos::View<double**>("V", n, m);
        auto T_view = Kokkos::View<double*>("T", lambda);
        auto result_view = Kokkos::View<double**>("result", n, m);
        
        // Create host mirrors and copy data
        auto V_host = Kokkos::create_mirror_view(V_view);
        auto T_host = Kokkos::create_mirror_view(T_view);
        
        // Copy input data
        for (int k = 0; k < m; k++) {
            for (int i = 0; i < n; i++) {
                V_host(i, k) = (*V)[i + k * n];
            }
        }
        for (int i = 0; i < lambda; i++) {
            T_host(i) = T[i];
        }
        
        // Copy to device
        Kokkos::deep_copy(V_view, V_host);
        Kokkos::deep_copy(T_view, T_host);
        
        // Main computation kernel
        Kokkos::parallel_for("ToeplitzMV", 
            Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, offset_edges}, {m, n - offset_edges}),
            KOKKOS_LAMBDA(const int k, const int i) {
                double sum = 0.0;
                const int j_first = (i - distcorrmin) > 0 ? (i - distcorrmin) : 0;
                const int j_last = (i + lambda) < n ? (i + lambda) : n;
                
                for (int j = j_first; j < j_last; j++) {
                    const int Tid = (j > i) ? (j - i) : (i - j);
                    if (Tid < lambda) {
                        sum += T_view(Tid) * V_view(j, k);
                    }
                }
                result_view(i, k) = sum;
            });
        
        // Copy result back to host
        auto result_host = Kokkos::create_mirror_view(result_view);
        Kokkos::deep_copy(result_host, result_view);
        
        // Copy back to original data structure
        for (int k = 0; k < m; k++) {
            for (int i = 0; i < n; i++) {
                (*V)[i + k * n] = result_host(i, k);
            }
        }
        
        return 0;
    } catch (...) {
        return -1;
    }
}

} // extern "C"