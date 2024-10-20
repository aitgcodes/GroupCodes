#!/bin/bash
#SBATCH --job-name=hf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 1200:00:00
##SBATCH --partition=LocalQ
#SBATCH --mem=20GB

cd $SLURM_SUBMIT_DIR

python3 hftraj.py >& log_hf_1
python3 hirsh.py  >& log_hf_2
bader density.cube
