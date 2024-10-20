#!/bin/bash
#SBATCH --job-name=hc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 1200:00:00
##SBATCH --partition=LocalQ
#SBATCH --mem=20GB

cd $SLURM_SUBMIT_DIR

DPATH=data

#t1='11.33e3 14.09e3 28.98e3 38.03e3'

export START_TIME=`date +%s.%3N`
export TIME_LIMIT=`squeue -h -o %L -j $SLURM_JOB_ID`

mpirun -np 1 --oversubscribe python3 src/unocc.py $DPATH/{gs/gs.gpw,unocc/{unocc{.out,.gpw},eig.dat,ksd{.out,.ulm}}}

mpirun -np 1 --oversubscribe python3 src/tdspec.py $DPATH/td-z-RPA/{dm.dat,abs.dat}

mpirun -np 1 --oversubscribe python3 src/pulse.py $DPATH/td-z-RPA/pulse/pulse.pickle 4.21

mpirun -np 8 --oversubscribe python3 src/frho.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/{fdm.ulm,frho}}

mpirun -np 8 --oversubscribe python3 src/pulse_dm.py $DPATH/td-z-RPA/{dm.dat,pulse/{pulse.pickle,dm.dat}} --maxtime 120e3 --freq_step 0.02

mpirun -np 8 --oversubscribe python3 src/pulse_rho.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/{frho,pulse/{pulse.pickle,rho}}} --maxtime 120e3 --freq_step 0.02

mpirun -np 8 --oversubscribe python3 src/pulse_energy.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/{pulse.pickle,rho,energy_all.npz}} --maxtime 120e3 --transitions all
mpirun -np 8 --oversubscribe python3 src/pulse_energy.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/{pulse.pickle,rho,energy_resonant.npz}} --maxtime 120e3 --transitions resonant
mpirun -np 8 --oversubscribe python3 src/pulse_energy.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/{pulse.pickle,rho,energy_nonresonant.npz}} --maxtime 120e3 --transitions nonresonant


