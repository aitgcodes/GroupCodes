#!/bin/bash
#SBATCH --job-name=hc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -t 1200:00:00
##SBATCH --partition=LocalQ
#SBATCH --mem=20GB

cd $SLURM_SUBMIT_DIR

DPATH=data

#time='11.4e3 11.42e3 11.44e3 11.46e3 11.5e3 11.54e3 17.3e3 17.32e3 17.34e3 17.4e3 17.44e3 25e3'
time='6e3 7e3 8e3 9e3 10e3 11e3'

export START_TIME=`date +%s.%3N`
export TIME_LIMIT=`squeue -h -o %L -j $SLURM_JOB_ID`

#mpirun -np 1 --oversubscribe python3 src/unocc.py $DPATH/{gs/gs.gpw,unocc/{unocc{.out,.gpw},eig.dat,ksd{.out,.ulm}}}
#
#mpirun -np 1 --oversubscribe python3 src/tdspec.py $DPATH/td-z-RPA/{dm.dat,abs.dat}
#
#mpirun -np 1 --oversubscribe python3 src/pulse.py $DPATH/td-z-RPA/pulse/pulse.pickle 4.14
#
#mpirun -np 8 --oversubscribe python3 src/frho.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/{fdm.ulm,frho}}
#
#mpirun -np 8 --oversubscribe python3 src/pulse_dm.py $DPATH/td-z-RPA/{dm.dat,pulse/{pulse.pickle,dm.dat}} --maxtime 30e3
#
#mpirun -np 8 --oversubscribe python3 src/pulse_rho.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/{frho,pulse/{pulse.pickle,rho}}} --maxtime 30e3
#
#mpirun -np 8 --oversubscribe python3 src/pulse_energy.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/{pulse.pickle,rho,energy_all.npz}} --maxtime 30e3 --transitions all
#mpirun -np 8 --oversubscribe python3 src/pulse_energy.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/{pulse.pickle,rho,energy_resonant.npz}} --maxtime 30e3 --transitions resonant
#mpirun -np 8 --oversubscribe python3 src/pulse_energy.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/{pulse.pickle,rho,energy_nonresonant.npz}} --maxtime 30e3 --transitions nonresonant

## Energy distributions as transition contribution maps
for t in $time
do
        mpirun -np 1 --oversubscribe python3 src/pulse_tcm.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/{pulse.pickle,rho,tcm_energy}} energy --time $t
done
#
## Hot-carrier distributions
for t in $time
do
        mpirun -np 1 --oversubscribe python3 src/pulse_hcdist.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/{pulse.pickle,rho,hcdist_all_$t.npz}}  --time $t  --transitions all
        mpirun -np 1 --oversubscribe python3 src/pulse_hcdist.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/{pulse.pickle,rho,hcdist_resonant_$t.npz}}  --time $t  --transitions resonant
done
#
#### Spatial atomic weights
#mpirun -np 16 --oversubscribe python3 src/atomweight.py $DPATH/unocc/{unocc.gpw,atomweight.npz} --eig_limits -6 6
#
#### Spatial hot-carrier density

for t in $time
do
        mpirun -np 16 --oversubscribe python3 src/pulse_hcdns.py $DPATH/{unocc/{unocc.gpw,ksd.ulm,atomweight.npz},td-z-RPA/pulse/{pulse.pickle,rho,he1dns_engt1eV_$t.cube}} electrons --time $t --energy_limits 1 10
        mpirun -np 16 --oversubscribe python3 src/pulse_hcdns.py $DPATH/{unocc/{unocc.gpw,ksd.ulm,atomweight.npz},td-z-RPA/pulse/{pulse.pickle,rho,hh1dns_$t.cube}} holes --time $t --energy_limits -10 0
done

# Generate cube files
#mpirun -np 1 --oversubscribe python3 src/generate_cubefiles.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/pulse.pickle} $DPATH/unocc/cube_files_resonant/ --transitions resonant 'gs'
#mpirun -np 1 --oversubscribe python3 src/generate_cubefiles.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/pulse.pickle} $DPATH/unocc/cube_files_nonresonant/ --transitions nonresonant 'gs'
##
## Generate weights of occupied and unoccupied Kohn-Sham states
#mpirun -np 1 --oversubscribe python3 src/weights_wfs.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/pulse.pickle} $DPATH/unocc/cube_files_resonant/ 'gs' --transitions resonant $DPATH/unocc/weights_resonant/{occ_np.dat,unocc_np.dat} --lx '-' --ly '-' --lz '31.7'
#mpirun -np 1 --oversubscribe python3 src/weights_wfs.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/pulse.pickle} $DPATH/unocc/cube_files_resonant/ 'gs' --transitions resonant $DPATH/unocc/weights_resonant/{occ_full.dat,unocc_full.dat} --lx '-' --ly '-' --lz '-'
#
#mpirun -np 1 --oversubscribe python3 src/weights_wfs.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/pulse.pickle} $DPATH/unocc/cube_files_nonresonant/ 'gs' --transitions nonresonant $DPATH/unocc/weights_nonresonant/{occ_np.dat,unocc_np.dat} --lx '-' --ly '-' --lz '31.7'
#mpirun -np 1 --oversubscribe python3 src/weights_wfs.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/pulse.pickle} $DPATH/unocc/cube_files_nonresonant/ 'gs' --transitions nonresonant $DPATH/unocc/weights_nonresonant/{occ_full.dat,unocc_full.dat} --lx '-' --ly '-' --lz '-'
##
#### Weighted transition probability
#mpirun -np 1 --oversubscribe python3 src/weighted_transp.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/{pulse.pickle,rho}} $DPATH/unocc/weights_resonant/{occ_full.dat,unocc_full.dat,occ_np.dat,unocc_np.dat}  $DPATH/td-z-RPA/pulse/wtp_resonant.npz  --maxtime 30e3  --transitions resonant
##
### Weighted transition probability distribution of partial processes
#mpirun -np 1 --oversubscribe python3 src/weighted_transp_dist.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/{pulse.pickle,rho}} $DPATH/unocc/weights_resonant/{occ_full.dat,unocc_full.dat,occ_np.dat,unocc_np.dat} $DPATH/td-z-RPA/pulse/wtp_dist_resonant_25.npz  --time 25e3  --transitions resonant
##
#### Weighted energy 
#mpirun -np 1 --oversubscribe python3 src/weighted_energy.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/{pulse.pickle,rho}} $DPATH/unocc/weights_nonresonant/{occ_full.dat,unocc_full.dat,occ_np.dat,unocc_np.dat}  $DPATH/td-z-RPA/pulse/wtenergy_nonresonant.npz  --maxtime 30e3  --transitions nonresonant

# Induced charge density
for t in $time
do
        mpirun -np 1 --oversubscribe python3 src/ind.py $DPATH/{unocc/{unocc.gpw,ksd.ulm},td-z-RPA/pulse/{rho,indden_t}} --time $t
done
