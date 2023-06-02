midapack_WF_env_alpha
# export DIR=/global/homes/m/mag/midapack/libmidapack


# cc test_spherical_harmonics_2.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HAT_DIR}/include -lcfitsio -L${S2HAT_DIR}/lib/cori/intel/ -ls2hat_std  -o test_functions/test_spherical_harmonics_2

# cc test_spherical_harmonics_2.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HAT_DIR}/include -lcfitsio -L${S2HAT_DIR}/lib/cori/intel/ -ls2hat_std   -L${HEALPIX_DIR}/lib/ -I${HEALPIX_DIR}/include/ -DHEALPIXDATA=${HEALPIX_DIR}share/healpix/  -o test_functions/test_spherical_harmonics_2

# cc test_functions/test_spherical_harmonics_2.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HAT_DIR}/include -lcfitsio -L${S2HAT_DIR}/lib/cori/intel/ -ls2hat_std   -L${HEALPIX_DIR}/lib/ -I${HEALPIX_DIR}/include/ -DHEALPIXDATA=${HEALPIX_DIR}share/healpix/  -o test_functions/test_spherical_harmonics_2
# cc test_functions/test_spherical_harmonics_2.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HAT_DIR}/include -lcfitsio -L${S2HAT_DIR}/lib/cori/intel/ -ls2hat_std   -L${HEALPIX_DIR}/lib/ -I${HEALPIX_DIR}/include/ -DHEALPIXDATA=${HEALPIX_DIR}share/healpix/  -o test_functions/test_spherical_harmonics_2b


# srun -n 2 test_functions/test_spherical_harmonics_2
# srun -n 1 test_functions/test_spherical_harmonics_2
# srun -n 1 test_functions/test_spherical_harmonics_2b

# cc test_wiener_filter/test_domain_generalization.c -Wall -I${PREFIX}/include/midapack -L${PREFIX}/lib -lmidapack -I${S2HAT_DIR}/include -lcfitsio -L${S2HAT_DIR}/lib/cori/intel/ -ls2hat_std   -L${HEALPIX_DIR}/lib/ -I${HEALPIX_DIR}/include/ -DHEALPIXDATA=${HEALPIX_DIR}share/healpix/  -o test_wiener_filter/test_domain_generalization
cc test_wiener_filter/test_domain_generalization.c -Wall -L${HEALPIX_DIR}/lib/ -I${HEALPIX_DIR}/include -DHEALPIXDATA=${HEALPIX_DIR}share/healpix/ -I${S2HAT_DIR}/include -lcfitsio -L${S2HAT_DIR}/lib/cori/intel/ -ls2hat_std -I${PREFIX}/include/midapack -I${PREFIX}/include/mappraiser -L${PREFIX}/lib -lmidapack -lmappraiser  -o test_wiener_filter/test_domain_generalization

srun -n 2 test_wiener_filter/test_domain_generalization
# srun -n 1 test_wiener_filter/test_domain_generalization // Need at least 2 for butterfly !


cc test_wiener_filter/map_to_white_noise.c -Wall -L${HEALPIX_DIR}/lib/ -I${HEALPIX_DIR}/include -DHEALPIXDATA=${HEALPIX_DIR}share/healpix/ -I${S2HAT_DIR}/include -lcfitsio -L${S2HAT_DIR}/lib/ -ls2hat_std -I${PREFIX}/include/midapack -I${PREFIX}/include/mappraiser -L${PREFIX}/lib -lmidapack -lmappraiser  -o test_wiener_filter/map_to_white_noise
ftn test_wiener_filter/map_to_white_noise.c -Wall -L${HEALPIX_DIR}/lib/ -I${HEALPIX_DIR}/include -DHEALPIXDATA=${HEALPIX_DIR}share/healpix/ -I${S2HAT_DIR}/include -L${S2HAT_DIR}/lib/ -ls2hat_std -I${PREFIX}/include/midapack -I${PREFIX}/include/mappraiser -L${PREFIX}/lib -lmidapack -lmappraiser  -o test_wiener_filter/map_to_white_noise


ftn test_wiener_filter/map_to_white_noise.c -Wall  -L${HEALPIX_DIR}/lib -I${HEALPIX_DIR}/include -I${S2HAT_DIR}/include -L${S2HAT_DIR}/lib/ -ls2hat_std -I${PREFIX}/include/midapack -I${PREFIX}/include/mappraiser -L${PREFIX}/lib -lmidapack -lmappraiser  -o test_wiener_filter/map_to_white_noise
srun -n 2 test_wiener_filter/map_to_white_noise

srun -n 2 test_wiener_filter/map_to_white_noise_0
srun -n 1 test_wiener_filter/map_to_white_noise_0
