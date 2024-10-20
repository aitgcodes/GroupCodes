#!/bin/bash
#SBATCH --job-name=ldos
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 1200:00:00
##SBATCH --partition=LocalQ
#SBATCH --mem=20GB

cd $SLURM_SUBMIT_DIR

/home/pramod/anaconda3/envs/gpaw_24.6.0/bin/mpirun -np 1 --oversubscribe /home/pramod/anaconda3/envs/gpaw_24.6.0/bin/python3 ldos.py >& ldos.log
/home/pramod/anaconda3/envs/gpaw_24.6.0/bin/mpirun -np 1 --oversubscribe /home/pramod/anaconda3/envs/gpaw_24.6.0/bin/python3 proj_dos_ao.py >& proj_dos_ao.log
