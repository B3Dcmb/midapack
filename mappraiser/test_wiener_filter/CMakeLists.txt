# Executable for running mappraiser without TOAST workflow

add_executable(map_to_white_noise_clean map_to_white_noise_clean.c)

set_property(TARGET map_to_white_noise_clean PROPERTY LINKER_LANGUAGE Fortran)

target_sources(map_to_white_noise_clean PRIVATE "map_to_white_noise_clean.c")
target_include_directories(map_to_white_noise_clean PRIVATE
        "../include"
        # "../include/mappraiser"
        "${CMAKE_SOURCE_DIR}/libmidapack/include"
        )

add_compile_options(${MPI_Fortran_COMPILE_FLAGS})

# MPI compile flag for midapack algebra modules
# target_compile_definitions(map_to_white_noise_clean PRIVATE W_MPI)

# if (HEALPIX_FOUND)
    #target_compile_definitions(s2hat_std PRIVATE HAVE_HEALPIX=1)
target_include_directories(map_to_white_noise_clean PUBLIC "${HEALPIX_INCLUDES}")
target_link_libraries(map_to_white_noise_clean PUBLIC "${HEALPIX_LIBRARIES}")
target_compile_options(map_to_white_noise_clean PUBLIC "${HEALPIX_LFLAGS}")
# endif (HEALPIX_FOUND)

# if (S2HAT_FOUND)
target_include_directories(map_to_white_noise_clean PRIVATE "${S2HAT_INCLUDES}")
target_link_libraries(map_to_white_noise_clean PRIVATE "${S2HAT_LIBRARIES}")
target_compile_options(map_to_white_noise_clean PRIVATE "${S2HAT_LFLAGS}")
# endif (S2HAT_FOUND)

# libraries
target_link_libraries(map_to_white_noise_clean PRIVATE mappraiser midapack MPI::MPI_C)
target_link_libraries(map_to_white_noise_clean PRIVATE MPI::MPI_Fortran)


install(TARGETS map_to_white_noise_clean 
        DESTINATION ${PROJECT_SOURCE_DIR}/mappraiser/test_wiener_filter/)
