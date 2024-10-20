#!/bin/bash

DPATH1=data
DPATH2=plots

at=/home/pramod/codes/hc/analysis-tools

t1='12.42e3 15.94e3 22.60e3 43.66e3 48.30e3 55.20e3 97.30e3 115.60e3'

# Energy as a function of time
#----------------RPA-----------
# Energy as a function of time
python3 $at/energy.py $DPATH1/td-z-RPA/pulse/energy_all.npz $DPATH2/energy/energy_all.dat
python3 $at/energy.py $DPATH1/td-z-RPA/pulse/energy_resonant.npz $DPATH2/energy/energy_resonant.dat
python3 $at/energy.py $DPATH1/td-z-RPA/pulse/energy_nonresonant.npz $DPATH2/energy/energy_nonresonant.dat
