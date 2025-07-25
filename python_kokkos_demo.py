#!/usr/bin/env python3
"""
Simple Python wrapper demonstrating Kokkos-accelerated MIDAPACK functionality.

This is a proof-of-concept for integrating Kokkos acceleration with Python.
For production use, this would be integrated with the existing Python bindings.
"""

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import os

class KokkosMidapack:
    """Python wrapper for Kokkos-accelerated MIDAPACK functions."""
    
    def __init__(self, library_path=None):
        """Initialize the wrapper and load the MIDAPACK library."""
        if library_path is None:
            # Try to find the library in common locations
            possible_paths = [
                './build/algebra/src/libmidapack.so',
                '/usr/local/lib/libmidapack.so',
                '/usr/lib/libmidapack.so'
            ]
            for path in possible_paths:
                if os.path.exists(path):
                    library_path = path
                    break
            
        if library_path is None:
            raise RuntimeError("Could not find libmidapack.so")
            
        self.lib = ctypes.CDLL(library_path)
        self._setup_function_signatures()
        self.kokkos_initialized = False
        
    def _setup_function_signatures(self):
        """Set up the function signatures for the C library."""
        # Kokkos initialization functions
        self.lib.kokkos_initialize.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)]
        self.lib.kokkos_initialize.restype = ctypes.c_int
        
        self.lib.kokkos_finalize.argtypes = []
        self.lib.kokkos_finalize.restype = None
        
        self.lib.kokkos_is_initialized.argtypes = []
        self.lib.kokkos_is_initialized.restype = ctypes.c_int
        
        # STMM function
        self.lib.stmm_simple_basic.argtypes = [
            ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),  # V
            ctypes.c_int,                                      # n
            ctypes.c_int,                                      # m  
            ctypes.POINTER(ctypes.c_double),                   # T
            ctypes.c_int,                                      # lambda
            ctypes.POINTER(ctypes.POINTER(ctypes.c_double))    # TV
        ]
        self.lib.stmm_simple_basic.restype = ctypes.c_int
        
    def initialize_kokkos(self):
        """Initialize Kokkos for GPU acceleration."""
        if self.kokkos_initialized:
            return
            
        # Create dummy command line arguments
        argc = 1
        argv = (ctypes.c_char_p * argc)(b"python")
        
        result = self.lib.kokkos_initialize(argc, argv)
        if result != 0:
            raise RuntimeError(f"Failed to initialize Kokkos (error code: {result})")
            
        self.kokkos_initialized = True
        print("Kokkos initialized successfully")
        
    def finalize_kokkos(self):
        """Finalize Kokkos."""
        if self.kokkos_initialized:
            self.lib.kokkos_finalize()
            self.kokkos_initialized = False
            print("Kokkos finalized")
            
    def is_kokkos_available(self):
        """Check if Kokkos is available and initialized."""
        return bool(self.lib.kokkos_is_initialized())
        
    def toeplitz_matvec(self, V, T):
        """
        Perform Toeplitz matrix-vector multiplication using Kokkos acceleration.
        
        Parameters:
        -----------
        V : numpy.ndarray
            Input matrix of shape (n, m)
        T : numpy.ndarray  
            Toeplitz coefficients of length lambda
            
        Returns:
        --------
        numpy.ndarray
            Result matrix of same shape as V
        """
        if not isinstance(V, np.ndarray) or not isinstance(T, np.ndarray):
            raise TypeError("V and T must be numpy arrays")
            
        if V.dtype != np.float64 or T.dtype != np.float64:
            raise TypeError("V and T must be float64 arrays")
            
        n, m = V.shape
        lambda_val = len(T)
        
        # Ensure arrays are C-contiguous
        V = np.ascontiguousarray(V.T)  # Transpose for C row-major layout
        T = np.ascontiguousarray(T)
        
        # Allocate output array
        TV = np.zeros_like(V)
        
        # Create pointers for the C function
        V_flat = V.flatten()
        TV_flat = TV.flatten()
        
        V_ptr = V_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        TV_ptr = TV_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        T_ptr = T.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        
        # Create double pointers
        V_double_ptr = ctypes.pointer(V_ptr)
        TV_double_ptr = ctypes.pointer(TV_ptr)
        
        # Call the C function
        result = self.lib.stmm_simple_basic(
            V_double_ptr, n, m, T_ptr, lambda_val, TV_double_ptr
        )
        
        if result != 0:
            raise RuntimeError(f"STMM computation failed (error code: {result})")
            
        # Reshape and transpose back to original layout
        return TV.reshape(m, n).T
        
    def __enter__(self):
        """Context manager entry."""
        self.initialize_kokkos()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.finalize_kokkos()


def demo_kokkos_acceleration():
    """Demonstration of Kokkos-accelerated computation from Python."""
    print("=== MIDAPACK Kokkos Python Demo ===\n")
    
    # Problem parameters
    n = 100      # Matrix size
    m = 5        # Number of columns  
    lambda_val = 7   # Bandwidth
    
    print(f"Problem size: n={n}, m={m}, lambda={lambda_val}")
    
    # Create test data
    V = np.random.randn(n, m)
    T = np.exp(-0.5 * np.arange(lambda_val))  # Exponential decay
    
    print(f"Input matrix V shape: {V.shape}")
    print(f"Toeplitz coefficients T: {T}")
    
    try:
        with KokkosMidapack() as midapack:
            print(f"Kokkos available: {midapack.is_kokkos_available()}")
            
            # Perform computation
            print("\nPerforming Toeplitz matrix-vector multiplication...")
            result = midapack.toeplitz_matvec(V, T)
            
            print(f"Result shape: {result.shape}")
            print(f"Result statistics:")
            print(f"  Min: {np.min(result):.6f}")
            print(f"  Max: {np.max(result):.6f}")
            print(f"  Mean: {np.mean(result):.6f}")
            print(f"  Std: {np.std(result):.6f}")
            
            print(f"\nFirst row of result: {result[0, :]}")
            
        print("\nDemo completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
        
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(demo_kokkos_acceleration())