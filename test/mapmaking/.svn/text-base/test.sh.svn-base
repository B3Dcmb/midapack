#!/bin/bash
#PBS -S /bin/bash
#PBS -N test1
#PBS -o out_test1
#PBS -e err_test1
#PBS -j oe
#PBS -m abe
#PBS -M fdauverg@apc.univ-paris7.fr
#PBS -l nodes=4:ppn=8,walltime=00:30:00
export SCRATCH="/scratch/$USER.$PBS_JOBID" # this set a system variable to point the temporary scratch repository
export PATH=/usr/local/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/openmpi/lib/:/usr/local/openmpi/lib/openmpi/:/usr/local/intel/fce/10.1.008/lib/:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
cd /home/dauvergn/midas-svn-OK/midas/trunk/test/mapmaking
mpirun -n 32 ./pcg0

