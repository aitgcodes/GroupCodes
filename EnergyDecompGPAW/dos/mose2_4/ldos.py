from gpaw import GPAW
import numpy as np
import os

directory='data'
if not os.path.exists(directory):
    os.makedirs(directory)

calc = GPAW('gs.gpw', txt=None)

ef = calc.get_fermi_level()

# Calculate the s, p, d projected LDOS for all the atoms

angular=['s', 'p', 'd']

# MoSe2 (N=4)
for a in range(46):
    for c in angular:
        energies, ldos = calc.get_orbital_ldos(a=a, spin=0, angular=c, npts=2001, width=0.1)
        outfile = open("data/tmd_{}_{}.dat".format(a, c), "w")
        for i in range(len(energies)):
            outfile.write(str(energies[i] - ef)+"    "+str(ldos[i])+"\n")
