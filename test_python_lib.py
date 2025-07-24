#!/usr/bin/env python3
"""
Simple Python test to verify that the Kokkos wrapper compiles and exports correctly.
"""

import ctypes
import os

def test_library_availability():
    """Test if the MIDAPACK library can be loaded and has Kokkos functions."""
    print("=== Testing MIDAPACK Kokkos Library Availability ===\n")
    
    # Try to find and load the library
    library_path = './build/algebra/src/libmidapack.so'
    
    if not os.path.exists(library_path):
        print(f"ERROR: Library not found at {library_path}")
        return 1
        
    try:
        lib = ctypes.CDLL(library_path)
        print(f"SUCCESS: Loaded library from {library_path}")
    except Exception as e:
        print(f"ERROR: Could not load library: {e}")
        return 1
    
    # Test if Kokkos functions are available
    kokkos_functions = [
        'kokkos_initialize',
        'kokkos_finalize', 
        'kokkos_is_initialized',
        'kokkos_stmm_simple_basic'
    ]
    
    available_functions = []
    missing_functions = []
    
    for func_name in kokkos_functions:
        try:
            func = getattr(lib, func_name)
            available_functions.append(func_name)
            print(f"✓ Found function: {func_name}")
        except AttributeError:
            missing_functions.append(func_name)
            print(f"✗ Missing function: {func_name}")
    
    print(f"\nSummary:")
    print(f"  Available functions: {len(available_functions)}")
    print(f"  Missing functions: {len(missing_functions)}")
    
    if missing_functions:
        print(f"  Missing: {missing_functions}")
        return 1
    else:
        print("  All Kokkos functions are available!")
        return 0

if __name__ == "__main__":
    import sys
    sys.exit(test_library_availability())