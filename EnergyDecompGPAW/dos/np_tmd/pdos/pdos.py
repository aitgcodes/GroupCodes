import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from ase.units import Hartree

from gpaw import GPAW

calc = GPAW('gs.gpw', txt=None)

ef = calc.get_fermi_level()

# Calculate the s, p, d projected LDOS for all the atoms

# Full DOS
list_0 = []
for a in range(101):
    list_0.append(a)

e, dos = calc.get_lcao_dos(atom_indices=list_0, basis_indices=None, npts=2001, width=0.1)

#e, dos = calc.get_dos(spin=0, npts=2001, width=0.1)

outfile_0 = open("Ag55_mose2_4_pdos.dat", "w")
for i in range(len(e)):
    outfile_0.write(str(e[i]-ef)+"    "+str(dos[i])+"\n")

# Ag55
list_1 = []
for a in range(55):
    list_1.append(a)

e, pdos = calc.get_lcao_dos(atom_indices=list_1, basis_indices=None, npts=2001, width=0.1)

outfile_1 = open("Ag55_pdos.dat", "w")
for i in range(len(e)):
    outfile_1.write(str(e[i]-ef)+"    "+str(pdos[i])+"\n")

# TMD
list_2 = []
for a in range(55, 101, 1):
    list_2.append(a)

e, pdos = calc.get_lcao_dos(atom_indices=list_2, basis_indices=None, npts=2001, width=0.1)
outfile_2 = open("MoSe2_4_pdos.dat", "w")
for i in range(len(e)):
    outfile_2.write(str(e[i]-ef)+"    "+str(pdos[i])+"\n")
