#!/bin/bash
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH --mail-user=magdy.morshed.fr@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 00:05:00

#OpenMP settings:
#export OMP_NUM_THREADS=1
#export OMP_PLACES=threads
#export OMP_PROC_BIND=spread

#virtual environnement
# conda activate CAMB-BB-py37

midapack_WF_env_alpha

cd /global/homes/m/mag/midapack/mappraiser/src/

cc test_wiener_filter/map_to_white_noise.c -Wall -L${HEALPIXROOT}/lib/ -I${HEALPIXROOT}/include -DHEALPIXDATA=${HEALPIXROOT}share/healpix/ -I${S2HATROOT}/include -lcfitsio -L${S2HATROOT}/lib/cori/intel/ -ls2hat_std -I${PREFIX}/include/midapack -I${PREFIX}/include/mappraiser -L${PREFIX}/lib -lmidapack -lmappraiser  -o test_wiener_filter/map_to_white_noise
srun -n 1 test_wiener_filter/map_to_white_noise


echo "Done !"
