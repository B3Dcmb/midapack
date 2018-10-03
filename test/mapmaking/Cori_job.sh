#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 2
#SBATCH -t 00:30:00
#SBATCH -J midapack_pcg
#SBATCH -C haswell
#SBATCH -L SCRATCH

source $HOME/.bashrc.ext
cd $SLURM_SUBMIT_DIR

time srun -n 128 ./pcg_pol > run.log 2>error.log
