#!/bin/sh
# -N 1       
#SBATCH --ntasks-per-node=48
#SBATCH --time=72:00:00
#SBATCH --job-name=np
#SBATCH --error=job.err
#SBATCH --output=job.out
#SBATCH --partition=standard

cd $SLURM_SUBMIT_DIR
DPATH=data

export START_TIME=`date +%s.%3N`
export TIME_LIMIT=`squeue -h -o %L -j $SLURM_JOB_ID`

export OMP_NUM_THREADS=1

time /home/pramodv.iitm/codes/GPAW/gpaw_24.6.0/bin/mpirun  -n $SLURM_NTASKS /home/pramodv.iitm/codes/GPAW/gpaw_24.6.0/bin/python3 src/gs.py structure/rlx-ico-Ag55.xyz $DPATH/gs/gs{.out,.gpw,-nonconv.gpw}
time /home/pramodv.iitm/codes/GPAW/gpaw_24.6.0/bin/mpirun  -n $SLURM_NTASKS /home/pramodv.iitm/codes/GPAW/gpaw_24.6.0/bin/python3 src/td.py $DPATH/{gs/gs.gpw,td-z-RPA/{td{.out,.gpw},dm.dat,fdm.ulm}} --fxc RPA --totaltime 120e3  --kick z --fdm_freq_step 0.02
