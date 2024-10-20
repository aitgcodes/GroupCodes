from ase.io import read
from ase.io import write
from ase.units import Bohr
from gpaw.analyse.hirshfeld import HirshfeldPartitioning
from gpaw import restart

atoms, calc = restart('gs.gpw')

# write Hirshfeld charges out
hf = HirshfeldPartitioning(calc)
for atom, charge in zip(atoms, hf.get_charges()):
    atom.charge = charge
# atoms.write('Hirshfeld.traj') # XXX Trajectory writer needs a fix
atoms.copy().write('Hirshfeld.traj')

# create electron density cube file ready for bader
rho = atoms.calc.get_all_electron_density(gridrefinement=4)
write('density.cube', atoms, data=rho * Bohr**3)
