# SPHERICAL HARMONICS
# to modify

target_sources(midapack PRIVATE
               "alm_tools.c"
               "covariance_matrix_tools.c"
               "files_io.c"
               "mpi_tools.c"
               "pix_tools.c"
               "s2hat_init_parameters.c"
               )


# if (MKL_FOUND)
#     target_compile_definitions(midapack PRIVATE HAVE_MKL)
#     target_include_directories(midapack PUBLIC "${MKL_INCLUDE_DIRS}")
#     target_link_libraries(midapack PUBLIC "${MKL_LIBRARIES}")
# elseif (LAPACK_FOUND)
#     # BLAS/LAPACK
#     #set_target_properties(mappraiser PROPERTIES LINK_FLAGS
#     #                      "${LAPACK_LINKER_FLAGS}" #${BLAS_LINKER_FLAGS}"
#     #                      )
#     #target_link_libraries(mappraiser PUBLIC "${LAPACK_LIBRARIES}")
#     #target_link_libraries(mappraiser "${BLAS_LIBRARIES}")
#     target_link_libraries(midapack PUBLIC LAPACK::LAPACK)
# endif ()
