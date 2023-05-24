# Executable for running mappraiser without TOAST workflow

add_executable(debug_run debug_run.cpp utils.c)

target_include_directories(debug_run PRIVATE
        "../include"
        "../include/mappraiser"
        "${CMAKE_SOURCE_DIR}/algebra/include"
        )

# MPI compile flag for midapack algebra modules
target_compile_definitions(debug_run PRIVATE W_MPI)

# libraries
target_link_libraries(debug_run PRIVATE mappraiser MPI::MPI_CXX)
