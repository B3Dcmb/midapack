#!/bin/bash -l
#SBATCH -q debug
#SBATCH -N 32
#SBATCH -t 00:30:00
#SBATCH -J midapack_map
#SBATCH -C haswell
#SBATCH -L SCRATCH
#DW jobdw capacity=3TiB access_mode=striped type=scratch
#DW stage_in source=/global/cscratch1/sd/elbouha/data_TOAST/test4_clean destination=$DW_JOB_STRIPED/test4_clean type=directory

source $HOME/.bashrc.ext
cd $SLURM_SUBMIT_DIR

export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OMP_NUM_THREADS=1

# TRIAL=2
# when the number of tasks is not a divisor of 64 use --cpu_bind = cores

# time srun -n 1 -c 2 --cpu_bind=cores ./toast_pcg > run.log 2>error.log
# time srun -n 16384 -c 2 --cpu_bind=cores ./toast_pcg_cyc $DW_JOB_STRIPED/test4_clean/ > run_16384.log 2>error_16384.log
# time srun -n 8192 -c 4 --cpu_bind=cores ./toast_pcg_cyc $DW_JOB_STRIPED/test4_clean/ > run_8192.log 2>error_8192.log
# time srun -n 4096 -c 8 --cpu_bind=cores ./toast_pcg_cyc $DW_JOB_STRIPED/test4_clean/ > run_4096.log 2>error_4096.log
# time srun -n 2048 -c 16 --cpu_bind=cores ./toast_pcg_cyc $DW_JOB_STRIPED/test4_clean/ > run_2048.log 2>error_2048.log
time srun -n 1024 -c 2 --cpu_bind=cores ./toast_pcg $DW_JOB_STRIPED/test4_clean/ > run.log 2>error.log
# time srun -n 512 -c 64 --cpu_bind=cores ./toast_pcg_cyc $DW_JOB_STRIPED/test4_clean/ > run_512.log 2>error_512.log
# time srun -n 256 -c 64 --cpu_bind=cores ./toast_pcg $DW_JOB_STRIPED/test3_clean/ > run_256.log 2>error_256.log
# time srun -n 128 -c 64 --cpu_bind=cores ./toast_pcg $DW_JOB_STRIPED/test3_clean/ > run_128.log 2>error_128.log


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
