midapack_WF_env_alpha
export DIR=/global/homes/m/mag/midapack
export MIDAPACK_LIB=${PREFIX}/midapack/lib
# INC='-I${DIR}/include -I${S2HATROOT}/include  -I${CFITSIOROOT}/include'
# LIB='-L${MIDAPACK_LIB} -lmidapack -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std  -L${CFITSIOROOT}/lib -lcfitsio'

### main_test_make_binary_mask
cc test_spherical_harmonics.c files_io.c -I${S2HATROOT}/include -I${CFITSIOROOT}/include -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -L${CFITSIOROOT}/lib -lcfitsio -o test_functions/test_spherical_harmonics


### main_S2HAT_GLOBAL_parameters
cc test_spherical_harmonics.c s2hat_init_parameters.c files_io.c -I${CFITSIOROOT}/include  -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std  -o test_functions/test_spherical_harmonics

# cc test_spherical_harmonics.c files_io.c s2hat_init_parameters.c -I${S2HATROOT}/include -I${CFITSIOROOT}/include -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -L${CFITSIOROOT}/lib -lcfitsio -o test_functions/test_spherical_harmonics


### main_S2HAT_LOCAL_parameters
# cc test_spherical_harmonics.c s2hat_init_parameters.c files_io.c $INC $LIB -o test_functions/test_spherical_harmonics

cc test_spherical_harmonics.c files_io.c s2hat_init_parameters.c -I${S2HATROOT}/include -I${CFITSIOROOT}/include -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -L${CFITSIOROOT}/lib -lcfitsio -o test_functions/test_spherical_harmonics

srun -n 2 test_functions/test_spherical_harmonics


### main_alm_pix_tools
cc test_spherical_harmonics.c files_io.c s2hat_init_parameters.c alm_pix_tools.c -I${S2HATROOT}/include -I${CFITSIOROOT}/include -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -L${CFITSIOROOT}/lib -lcfitsio -o test_functions/test_spherical_harmonics

srun -n 1 test_functions/test_spherical_harmonics



### main_covariance_matrix_tools
cc -D W_MKL test_spherical_harmonics.c files_io.c s2hat_init_parameters.c alm_pix_tools.c covariance_matrix_tools.c -I${S2HATROOT}/include -I${CFITSIOROOT}/include -I${MKLROOT}/include  -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -L${CFITSIOROOT}/lib -lcfitsio -L${MKLROOT}/lib/intel64 -lmkl_rt -o test_functions/test_spherical_harmonics

srun -n 1 test_functions/test_spherical_harmonics

cc -D W_MKL -D W_MPI test_spherical_harmonics.c files_io.c s2hat_init_parameters.c alm_pix_tools.c covariance_matrix_tools.c -I${S2HATROOT}/include -I${CFITSIOROOT}/include -I${MKLROOT}/include -I${DIR}/include  -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -L${CFITSIOROOT}/lib -lcfitsio -L${MKLROOT}/lib/intel64 -lmkl_rt -L${MIDAPACK_LIB} -lmidapack -o test_functions/test_spherical_harmonics

srun -n 2 test_functions/test_spherical_harmonics




cc -D W_MKL -D W_MPI test_spherical_harmonics.c files_io.c s2hat_init_parameters.c alm_pix_tools.c covariance_matrix_tools.c tmp_communication.c -I${S2HATROOT}/include -I${CFITSIOROOT}/include -I${MKLROOT}/include -I${DIR}/include  -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -L${CFITSIOROOT}/lib -lcfitsio -L${MKLROOT}/lib/intel64 -lmkl_rt -L${MIDAPACK_LIB} -lmidapack -o test_functions/test_spherical_harmonics

# srun -n 1 test_functions/test_spherical_harmonics


#### OLD, not working anyway
# all_obj = $(WIENERFILTER_OBJ)/alm_pix_tools.o $(WIENERFILTER_OBJ)/covariance_matrix_tools.o $(WIENERFILTER_OBJ)/s2hat_init_parameters.o $(WIENERFILTER_OBJ)/files_io.o $(WIENERFILTER_OBJ)/tmp_communication.o
# test_spherical_harmonics.o: test_spherical_harmonics.c
# 	$(CC) -c $< $(INC) $(LIB) $(CFITSIO_INC) $(CFITSIO_LIB) -o $@

# test_spherical_harmonics: test_spherical_harmonics.o $(all_obj)
# 	$(CC) $(INC) $(LIB) $(CFITSIO_INC) $(CFITSIO_LIB) $(all_obj) -o $@
