midapack_WF_env_alpha
export DIR=/global/homes/m/mag/midapack/libmidapack
# export MIDAPACK_LIB=${PREFIX}/midapack/lib
# INC='-I${DIR}/include -I${S2HATROOT}/include  -I${CFITSIOROOT}/include'
# LIB='-L${MIDAPACK_LIB} -lmidapack -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std  -L${CFITSIOROOT}/lib -lcfitsio'
INC='-I${PREFIX}/include/midapack -libmidapack -I${S2HATROOT}/include'
LIB='-lmidapack -L${PREFIX}/lib -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std'
### main_test_make_binary_mask
# cc test_spherical_harmonics_2.c  -I${S2HATROOT}/include -I${CFITSIOROOT}/include -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -L${CFITSIOROOT}/lib -lcfitsio -o test_functions/test_spherical_harmonics_2

cc test_spherical_harmonics_2.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std  -o test_functions/test_spherical_harmonics_2





#### OLD
cc test_spherical_harmonics_2.c  -I${PREFIX}/include/midapack -lmidapack -L${PREFIX}/lib -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std  -o test_functions/test_spherical_harmonics_2


### main_S2HAT_GLOBAL_parameters
cc test_spherical_harmonics_2.c s2hat_init_parameters.c files_io.c -I${CFITSIOROOT}/include  -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std  -o test_functions/test_spherical_harmonics_2

# cc test_spherical_harmonics_2.c files_io.c s2hat_init_parameters.c -I${S2HATROOT}/include -I${CFITSIOROOT}/include -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -L${CFITSIOROOT}/lib -lcfitsio -o test_functions/test_spherical_harmonics_2


### main_S2HAT_LOCAL_parameters
# cc test_spherical_harmonics_2.c s2hat_init_parameters.c files_io.c $INC $LIB -o test_functions/test_spherical_harmonics_2

cc test_spherical_harmonics_2.c files_io.c s2hat_init_parameters.c -I${S2HATROOT}/include -I${CFITSIOROOT}/include -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -L${CFITSIOROOT}/lib -lcfitsio -o test_functions/test_spherical_harmonics_2

srun -n 2 test_functions/test_spherical_harmonics_2


### main_alm_pix_tools
cc test_spherical_harmonics_2.c files_io.c s2hat_init_parameters.c alm_pix_tools.c -I${S2HATROOT}/include -I${CFITSIOROOT}/include -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -L${CFITSIOROOT}/lib -lcfitsio -o test_functions/test_spherical_harmonics_2

srun -n 1 test_functions/test_spherical_harmonics_2



### main_covariance_matrix_tools
cc -D W_MKL test_spherical_harmonics_2.c files_io.c s2hat_init_parameters.c alm_pix_tools.c covariance_matrix_tools.c -I${S2HATROOT}/include -I${CFITSIOROOT}/include -I${MKLROOT}/include  -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -L${CFITSIOROOT}/lib -lcfitsio -L${MKLROOT}/lib/intel64 -lmkl_rt -o test_functions/test_spherical_harmonics_2

srun -n 1 test_functions/test_spherical_harmonics_2

cc -D W_MKL -D W_MPI test_spherical_harmonics_2.c files_io.c s2hat_init_parameters.c alm_pix_tools.c covariance_matrix_tools.c -I${S2HATROOT}/include -I${CFITSIOROOT}/include -I${MKLROOT}/include -I${DIR}/include  -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -L${CFITSIOROOT}/lib -lcfitsio -L${MKLROOT}/lib/intel64 -lmkl_rt -L${MIDAPACK_LIB} -lmidapack -o test_functions/test_spherical_harmonics_2

srun -n 2 test_functions/test_spherical_harmonics_2




cc -D W_MKL -D W_MPI test_spherical_harmonics_2.c files_io.c s2hat_init_parameters.c alm_pix_tools.c covariance_matrix_tools.c tmp_communication.c -I${S2HATROOT}/include -I${CFITSIOROOT}/include -I${MKLROOT}/include -I${DIR}/include  -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -L${CFITSIOROOT}/lib -lcfitsio -L${MKLROOT}/lib/intel64 -lmkl_rt -L${MIDAPACK_LIB} -lmidapack -o test_functions/test_spherical_harmonics_2

# srun -n 1 test_functions/test_spherical_harmonics_2


#### OLD, not working anyway
# all_obj = $(WIENERFILTER_OBJ)/alm_pix_tools.o $(WIENERFILTER_OBJ)/covariance_matrix_tools.o $(WIENERFILTER_OBJ)/s2hat_init_parameters.o $(WIENERFILTER_OBJ)/files_io.o $(WIENERFILTER_OBJ)/tmp_communication.o
# test_spherical_harmonics_2.o: test_spherical_harmonics_2.c
# 	$(CC) -c $< $(INC) $(LIB) $(CFITSIO_INC) $(CFITSIO_LIB) -o $@

# test_spherical_harmonics_2: test_spherical_harmonics_2.o $(all_obj)
# 	$(CC) $(INC) $(LIB) $(CFITSIO_INC) $(CFITSIO_LIB) $(all_obj) -o $@
