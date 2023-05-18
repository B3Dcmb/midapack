# ECG segmentation fault debug

add_executable(ecg_debug utils.c utils.h)

target_sources(ecg_debug PRIVATE "ecg_debug.c")
target_include_directories(ecg_debug PRIVATE
        "../../include"
        "../../include/mappraiser"
        "${CMAKE_SOURCE_DIR}/algebra/include"
        )

# MPI compile flag for midapack algebra modules
target_compile_definitions(ecg_debug PRIVATE W_MPI)

# libraries
target_link_libraries(ecg_debug PRIVATE mappraiser MPI::MPI_C)
