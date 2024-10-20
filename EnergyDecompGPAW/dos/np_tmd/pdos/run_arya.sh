#!/bin/bash
#SBATCH --job-name=pdos
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 1200:00:00
##SBATCH --partition=LocalQ
#SBATCH --mem=20GB

cd $SLURM_SUBMIT_DIR

/home/pramod/anaconda3/envs/gpaw_24.6.0/bin/mpirun -np 1 --oversubscribe /home/pramod/anaconda3/envs/gpaw_24.6.0/bin/python3 pdos.py >& pdos.log
