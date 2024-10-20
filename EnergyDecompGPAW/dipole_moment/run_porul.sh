#!/bin/sh
# -N 1       
#SBATCH --ntasks-per-node=48
#SBATCH --time=72:00:00
#SBATCH --job-name=np
#SBATCH --error=job.err
#SBATCH --output=job.out
#SBATCH --partition=standard

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=1

time /home/pramodv.iitm/codes/GPAW/gpaw_24.6.0/bin/mpirun  -n $SLURM_NTASKS /home/pramodv.iitm/codes/GPAW/gpaw_24.6.0/bin/python3 gs.py 
time /home/pramodv.iitm/codes/GPAW/gpaw_24.6.0/bin/mpirun  -n $SLURM_NTASKS /home/pramodv.iitm/codes/GPAW/gpaw_24.6.0/bin/python3 td.py 
