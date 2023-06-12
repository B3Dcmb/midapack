# Executable for running mappraiser without TOAST workflow

add_executable(mini_test_white_noise_1 mini_test_white_noise_1.c)

set_property(TARGET mini_test_white_noise_1 PROPERTY LINKER_LANGUAGE Fortran)

target_sources(mini_test_white_noise_1 PRIVATE "mini_test_white_noise_1.c")
target_include_directories(mini_test_white_noise_1 PRIVATE
        "../include"
        # "../include/mappraiser"
        "${CMAKE_SOURCE_DIR}/libmidapack/include"
        )

add_compile_options(${MPI_Fortran_COMPILE_FLAGS})

# MPI compile flag for midapack algebra modules
# target_compile_definitions(mini_test_white_noise_1 PRIVATE W_MPI)

# if (HEALPIX_FOUND)
    #target_compile_definitions(s2hat_std PRIVATE HAVE_HEALPIX=1)
target_include_directories(mini_test_white_noise_1 PUBLIC "${HEALPIX_INCLUDES}")
target_link_libraries(mini_test_white_noise_1 PUBLIC "${HEALPIX_LIBRARIES}")
target_compile_options(mini_test_white_noise_1 PUBLIC "${HEALPIX_LFLAGS}")
# endif (HEALPIX_FOUND)

# if (S2HAT_FOUND)
target_include_directories(mini_test_white_noise_1 PRIVATE "${S2HAT_INCLUDES}")
target_link_libraries(mini_test_white_noise_1 PRIVATE "${S2HAT_LIBRARIES}")
target_compile_options(mini_test_white_noise_1 PRIVATE "${S2HAT_LFLAGS}")
# endif (S2HAT_FOUND)


target_include_directories(mini_test_white_noise_1 PRIVATE "${CFITSIO_INCLUDES}")
target_link_libraries(mini_test_white_noise_1 PUBLIC "${CFITSIO_LIBRARIES}")

# libraries
target_link_libraries(mini_test_white_noise_1 PRIVATE mappraiser midapack MPI::MPI_C)
target_link_libraries(mini_test_white_noise_1 PRIVATE MPI::MPI_Fortran)


install(TARGETS mini_test_white_noise_1 
        DESTINATION ${PROJECT_SOURCE_DIR}/mappraiser/test_wiener_filter/)
