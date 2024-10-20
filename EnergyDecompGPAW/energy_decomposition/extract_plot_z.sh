#!/bin/bash

DPATH1=data
DPATH2=plots/td-z-RPA

at=/home/pramod/codes/hc/analysis-tools

time='6e3 7e3 8e3 9e3 10e3 11e3 11.4e3 11.42e3 11.44e3 11.46e3 11.5e3 11.54e3 17.3e3 17.32e3 17.34e3 17.4e3 17.44e3 25e3'
time_tcm='06000 07000 08000 09000 10000 11000 11400 11420 11440 11460 11500 11540 17300 17320 17340 17400 17440 25000'
### Energy decomposition
#python3 $at/energy.py $DPATH1/td-z-RPA/pulse/energy_all.npz $DPATH2/energy/energy_all.dat
#python3 $at/energy.py $DPATH1/td-z-RPA/pulse/energy_resonant.npz $DPATH2/energy/energy_resonant.dat
#python3 $at/energy.py $DPATH1/td-z-RPA/pulse/energy_nonresonant.npz $DPATH2/energy/energy_nonresonant.dat
##
#### Weighted transition energy
#python3 $at/wtenergy.py $DPATH1/td-z-RPA/pulse/wtenergy_nonresonant.npz $DPATH2/wtenergy/wtenergy_nonresonant.dat
##
#### Weighted transition probability
#python3 $at/wtprob.py $DPATH1/td-z-RPA/pulse/wtp_resonant.npz $DPATH2/wtp/wtp_resonant.dat $DPATH2/wtp/total_prob.dat
##
#### Weighted transition probability distribution
#python3 $at/wtpdist.py $DPATH1/td-z-RPA/pulse/wtp_dist_resonant_25.npz $DPATH2/wtp_dist/wtp_dist_resonant_25.dat
#
#---Hot-carrier distribution---------------------------------------------
for t in $time
do
   python3 $at/hcdist.py $DPATH1/td-z-RPA/pulse/hcdist_all_$t.npz  $DPATH2/hcdist/hcdist_all_$t.dat
   python3 $at/hcdist.py $DPATH1/td-z-RPA/pulse/hcdist_resonant_$t.npz  $DPATH2/hcdist/hcdist_resonant_$t.dat
done
#### Energy as transition contribution map
for t in $time_tcm
do
        python3 $at/plot_tcm.py $DPATH1/unocc/unocc.gpw $DPATH1/unocc/ksd.ulm  $DPATH1/td-z-RPA/pulse/tcm_energy/t00$t.0.npz $DPATH2/tcm/energy --freq 4.14
done
