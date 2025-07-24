/**
 * @file kokkos_wrapper.h
 * @brief C wrapper for Kokkos functionality
 * @author Midapack Team
 * @date 2024
 * 
 * This header provides a C interface to Kokkos functionality,
 * allowing the existing C codebase to use Kokkos for GPU/CPU acceleration.
 */

#ifndef KOKKOS_WRAPPER_H
#define KOKKOS_WRAPPER_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Initialize Kokkos
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return 0 on success, non-zero on error
 */
int kokkos_initialize(int argc, char** argv);

/**
 * @brief Finalize Kokkos
 */
void kokkos_finalize(void);

/**
 * @brief Check if Kokkos is initialized
 * @return 1 if initialized, 0 otherwise
 */
int kokkos_is_initialized(void);

/**
 * @brief Allocate memory on device
 * @param size Size in bytes
 * @return Pointer to device memory
 */
void* kokkos_malloc(size_t size);

/**
 * @brief Free device memory
 * @param ptr Pointer to device memory
 */
void kokkos_free(void* ptr);

/**
 * @brief Copy data from host to device
 * @param dst Device pointer
 * @param src Host pointer
 * @param size Size in bytes
 */
void kokkos_memcpy_to_device(void* dst, const void* src, size_t size);

/**
 * @brief Copy data from device to host
 * @param dst Host pointer
 * @param src Device pointer
 * @param size Size in bytes
 */
void kokkos_memcpy_to_host(void* dst, const void* src, size_t size);

/**
 * @brief Kokkos-accelerated Toeplitz matrix-vector multiplication
 * @param V Input/output data matrix
 * @param n Number of rows
 * @param m Number of columns
 * @param T Toeplitz matrix data
 * @param lambda Band width
 * @return 0 on success, non-zero on error
 */
int kokkos_stmm_simple_basic(double** V, int n, int m, const double* T, int lambda);

#ifdef __cplusplus
}
#endif

#endif /* KOKKOS_WRAPPER_H */