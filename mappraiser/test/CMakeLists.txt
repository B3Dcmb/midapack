# ________________________________________________________
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

# ____________________________________________________
# Test executable for baseline computation performance

add_executable(baseline baseline.cpp)

target_include_directories(baseline PRIVATE
        "../include"
        "../include/mappraiser"
)

# libraries
target_link_libraries(baseline PRIVATE mappraiser OpenMP::OpenMP_CXX)

# ____________________________________________________
# Test executable for mapmat testing

add_executable(mapmat mapmat.cpp)

target_include_directories(mapmat PRIVATE
        "../include"
        "../include/mappraiser"
        "../../algebra/include"
        "../../algebra/include/mapmat"
)

target_compile_definitions(mapmat PRIVATE W_MPI)
target_link_libraries(mapmat PRIVATE midapack MPI::MPI_CXX)

# ____________________________________________________
# Test executable for Gap testing

add_executable(gaps mapping_gaps.cpp)

target_include_directories(gaps PRIVATE
        "../include"
        "../include/mappraiser"
        "../../algebra/include"
        "../../algebra/include/mapmat"
)

target_compile_definitions(gaps PRIVATE W_MPI)
target_link_libraries(gaps PRIVATE mappraiser MPI::MPI_CXX)

# ____________________________________________________
# Test executable for shift-and-overlap testing

add_executable(shift_overlap shift_overlap.cpp utils.c)

target_include_directories(shift_overlap PRIVATE
        "../include"
        "../include/mappraiser"
        "../../algebra/include"
        "../../algebra/include/mapmat"
)

target_compile_definitions(shift_overlap PRIVATE W_MPI)
target_link_libraries(shift_overlap PRIVATE mappraiser MPI::MPI_CXX)
