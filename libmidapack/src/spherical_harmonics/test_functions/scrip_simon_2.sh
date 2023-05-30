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

cc test_spherical_harmonics_2.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std   -L${HEALPIXROOT}/lib/ -I${HEALPIXROOT}/include/ -DHEALPIXDATA=${HEALPIXROOT}share/healpix/  -o test_functions/test_spherical_harmonics_2

cc test_functions/test_spherical_harmonics_2.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std   -L${HEALPIXROOT}/lib/ -I${HEALPIXROOT}/include/ -DHEALPIXDATA=${HEALPIXROOT}share/healpix/  -o test_functions/test_spherical_harmonics_2
cc test_functions/test_spherical_harmonics_2.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std   -L${HEALPIXROOT}/lib/ -I${HEALPIXROOT}/include/ -DHEALPIXDATA=${HEALPIXROOT}share/healpix/  -o test_functions/test_spherical_harmonics_2b


srun -n 2 test_functions/test_spherical_harmonics_2
srun -n 1 test_functions/test_spherical_harmonics_2
srun -n 1 test_functions/test_spherical_harmonics_2b

cc test_functions/test_butterfly.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std   -L${HEALPIXROOT}/lib/ -I${HEALPIXROOT}/include/ -DHEALPIXDATA=${HEALPIXROOT}share/healpix/  -o test_functions/test_butterfly

srun -n 2 test_functions/test_butterfly


cc test_functions/test_spherical_harmonics_2b.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std   -L${HEALPIXROOT}/lib/ -I${HEALPIXROOT}/include/ -DHEALPIXDATA=${HEALPIXROOT}share/healpix/  -o test_functions/test_spherical_harmonics_2b

srun -n 2 test_functions/test_spherical_harmonics_2b
srun -n 1 test_functions/test_spherical_harmonics_2b


#### Test base funcs midapack ####

cc test_functions/test_butterfly_1.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std   -L${HEALPIXROOT}/lib/ -I${HEALPIXROOT}/include/ -DHEALPIXDATA=${HEALPIXROOT}share/healpix/  -o test_functions/test_butterfly_1

srun -n 2 test_functions/test_butterfly_1


cc test_functions/test_als_0.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std   -L${HEALPIXROOT}/lib/ -I${HEALPIXROOT}/include/ -DHEALPIXDATA=${HEALPIXROOT}share/healpix/  -o test_functions/test_als_0

srun -n 1 test_functions/test_als_0
