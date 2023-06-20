# Executable for running mappraiser without TOAST workflow

add_executable(cut_sky_test_white_noise_diag cut_sky_test_white_noise_diag.c)

set_property(TARGET cut_sky_test_white_noise_diag PROPERTY LINKER_LANGUAGE Fortran)

target_sources(cut_sky_test_white_noise_diag PRIVATE "cut_sky_test_white_noise_diag.c")
target_include_directories(cut_sky_test_white_noise_diag PRIVATE
        "../include"
        # "../include/mappraiser"
        "${CMAKE_SOURCE_DIR}/libmidapack/include"
        )

add_compile_options(${MPI_Fortran_COMPILE_FLAGS})

# MPI compile flag for midapack algebra modules
# target_compile_definitions(cut_sky_test_white_noise_diag PRIVATE W_MPI)

# if (HEALPIX_FOUND)
    #target_compile_definitions(s2hat_std PRIVATE HAVE_HEALPIX=1)
target_include_directories(cut_sky_test_white_noise_diag PUBLIC "${HEALPIX_INCLUDES}")
target_link_libraries(cut_sky_test_white_noise_diag PUBLIC "${HEALPIX_LIBRARIES}")
target_compile_options(cut_sky_test_white_noise_diag PUBLIC "${HEALPIX_LFLAGS}")
# endif (HEALPIX_FOUND)

# if (S2HAT_FOUND)
target_include_directories(cut_sky_test_white_noise_diag PRIVATE "${S2HAT_INCLUDES}")
target_link_libraries(cut_sky_test_white_noise_diag PRIVATE "${S2HAT_LIBRARIES}")
target_compile_options(cut_sky_test_white_noise_diag PRIVATE "${S2HAT_LFLAGS}")
# endif (S2HAT_FOUND)


target_include_directories(cut_sky_test_white_noise_diag PRIVATE "${CFITSIO_INCLUDES}")
target_link_libraries(cut_sky_test_white_noise_diag PUBLIC "${CFITSIO_LIBRARIES}")

# libraries
target_link_libraries(cut_sky_test_white_noise_diag PRIVATE mappraiser midapack MPI::MPI_C)
target_link_libraries(cut_sky_test_white_noise_diag PRIVATE MPI::MPI_Fortran)


install(TARGETS cut_sky_test_white_noise_diag 
        DESTINATION ${PROJECT_SOURCE_DIR}/mappraiser/test_wiener_filter/)
