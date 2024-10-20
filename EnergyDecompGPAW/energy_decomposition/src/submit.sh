#!/bin/bash
#SBATCH -A account

export START_TIME=`date +%s.%3N`
export TIME_LIMIT=`squeue -h -o %L -j $SLURM_JOB_ID`

module load gpaw/version
export GPAW_SETUP_PATH=gpaw-setups

echo "$*"
srun gpaw python $*
