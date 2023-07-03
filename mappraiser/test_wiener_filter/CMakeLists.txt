# Executable for running mappraiser without TOAST workflow
# include(CTest)
add_executable(test_WF test_WF.c)
# add_test(test_WF_check test_WF)

# find_program( MEMORYCHECK_COMMAND valgrind )
# set( MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full" )
# set( MEMORYCHECK_SUPPRESSIONS_FILE "${PROJECT_SOURCE_DIR}/valgrind_suppress.txt" )

set_property(TARGET test_WF PROPERTY LINKER_LANGUAGE Fortran)

target_sources(test_WF PRIVATE "test_WF.c")
target_include_directories(test_WF PRIVATE
        "../include"
        # "../include/mappraiser"
        "${CMAKE_SOURCE_DIR}/algebra/include"
        )

add_compile_options(${MPI_Fortran_COMPILE_FLAGS})

# MPI compile flag for midapack algebra modules
# target_compile_definitions(test_WF PRIVATE W_MPI)

# if (HEALPIX_FOUND)
    #target_compile_definitions(s2hat_std PRIVATE HAVE_HEALPIX=1)
target_include_directories(test_WF PUBLIC "${HEALPIX_INCLUDES}")
target_link_libraries(test_WF PUBLIC "${HEALPIX_LIBRARIES}")
target_compile_options(test_WF PUBLIC "${HEALPIX_LFLAGS}")
# endif (HEALPIX_FOUND)

# if (S2HAT_FOUND)
target_include_directories(test_WF PRIVATE "${S2HAT_INCLUDES}")
target_link_libraries(test_WF PRIVATE "${S2HAT_LIBRARIES}")
target_compile_options(test_WF PRIVATE "${S2HAT_LFLAGS}")
# endif (S2HAT_FOUND)


target_include_directories(test_WF PRIVATE "${CFITSIO_INCLUDES}")
target_link_libraries(test_WF PUBLIC "${CFITSIO_LIBRARIES}")

# libraries
target_link_libraries(test_WF PRIVATE mappraiser midapack MPI::MPI_C)
target_link_libraries(test_WF PRIVATE MPI::MPI_Fortran)


install(TARGETS test_WF 
        DESTINATION ${PROJECT_SOURCE_DIR}/mappraiser/test_wiener_filter/)
