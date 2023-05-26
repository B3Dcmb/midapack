midapack_WF_env_alpha
# export DIR=/global/homes/m/mag/midapack/libmidapack


# cc test_spherical_harmonics_2.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std  -o test_functions/test_spherical_harmonics_2

# cc test_spherical_harmonics_2.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std   -L${HEALPIXROOT}/lib/ -I${HEALPIXROOT}/include/ -DHEALPIXDATA=${HEALPIXROOT}share/healpix/  -o test_functions/test_spherical_harmonics_2

# cc test_functions/test_spherical_harmonics_2.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std   -L${HEALPIXROOT}/lib/ -I${HEALPIXROOT}/include/ -DHEALPIXDATA=${HEALPIXROOT}share/healpix/  -o test_functions/test_spherical_harmonics_2
# cc test_functions/test_spherical_harmonics_2.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std   -L${HEALPIXROOT}/lib/ -I${HEALPIXROOT}/include/ -DHEALPIXDATA=${HEALPIXROOT}share/healpix/  -o test_functions/test_spherical_harmonics_2b


# srun -n 2 test_functions/test_spherical_harmonics_2
# srun -n 1 test_functions/test_spherical_harmonics_2
# srun -n 1 test_functions/test_spherical_harmonics_2b

# cc test_wiener_filter/test_domain_generalization.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std   -L${HEALPIXROOT}/lib/ -I${HEALPIXROOT}/include/ -DHEALPIXDATA=${HEALPIXROOT}share/healpix/  -o test_wiener_filter/test_domain_generalization
cc test_wiener_filter/test_domain_generalization.c -Wall -L${HEALPIXROOT}/lib/ -I${HEALPIXROOT}/include -DHEALPIXDATA=${HEALPIXROOT}share/healpix/ -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -I${PREFIX}/include/midapack -I${PREFIX}/include/mappraiser -L${PREFIX}/lib -lmidapack -lmappraiser  -o test_wiener_filter/test_domain_generalization

srun -n 2 test_wiener_filter/test_domain_generalization
# srun -n 1 test_wiener_filter/test_domain_generalization // Need at least 2 for butterfly !
