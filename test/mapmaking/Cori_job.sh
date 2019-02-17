#!/bin/bash -l
#SBATCH -q debug
#SBATCH -N 1
#SBATCH -t 00:05:00
#SBATCH -J midapack_map
#SBATCH -C haswell
#SBATCH -L SCRATCH

source $HOME/.bashrc.ext
cd $SLURM_SUBMIT_DIR

export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NUM_THREADS=2

# TRIAL=2
# when the number of tasks is not a divisor of 64 use --cpu_bind = cores

time srun -n 16 -c 4 ./toast_pcg > run.log 2>error.log
# time srun -n 128 ./pcg_pol > run.log 2>error.log
# time srun -n 32 -N 512 -c 1 ./test_com > run_32_${TRIAL}.log 2>error_32_${TRIAL}.log &
# time srun -n 64 -N 512 -c 1 ./test_com > run_64_${TRIAL}.log 2>error_64_${TRIAL}.log &
# time srun -n 128 -N 512 -c 1 ./test_com > run_128_${TRIAL}.log 2>error_128_${TRIAL}.log &
# time srun -n 256 -N 512 -c 1 ./test_com > run_256_${TRIAL}.log 2>error_256_${TRIAL}.log &
# time srun -n 512 -N 512 -c 1 ./test_com > run_512_${TRIAL}.log 2>error_512_${TRIAL}.log &
# time srun -n 1024 -N 512 -c 1 ./test_com > run_1024_${TRIAL}.log 2>error_1024_${TRIAL}.log &
# time srun -n 2048 -N 512 -c 1 ./test_com > run_2048_${TRIAL}.log 2>error_2048_${TRIAL}.log &
# time srun -n 4096 -N 512 -c 1 ./test_com > run_4096_${TRIAL}.log 2>error_4096_${TRIAL}.log &
# time srun -n 8192 -N 512 -c 1 ./test_com > run_8192_${TRIAL}.log 2>error_8192_${TRIAL}.log &
# time srun -n 16384 -N 512 -c 1 ./test_com > run_16384_${TRIAL}.log 2>error_16384_${TRIAL}.log &
# wait
