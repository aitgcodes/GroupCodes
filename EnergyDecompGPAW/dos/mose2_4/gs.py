import os
import sys
from ase.io import read
from ase.parallel import paropen
from gpaw import GPAW, FermiDirac, PoissonSolver, Mixer
from gpaw import KohnShamConvergenceError
from gpaw.cluster import Cluster
from gpaw.utilities.timelimit import TimeLimiter, time_to_seconds
from gpaw.mpi import world
from parallel_util import get_parallel
from gpaw.poisson_moment import MomentCorrectionPoissonSolver

# Read the structure from the xyz file
atoms = read('rlx-mose2_4.xyz')
atoms.center(vacuum=6.0)

# Increase the accuracy of density for ground state
convergence = {'density': 1e-12}

# Use occupation smearing, weak mixer and GLLB weight smearing
# to facilitate convergence
occupations = FermiDirac(width=0.05)
mixer=Mixer(0.1, 5, 1)
xc = 'GLLBSC'
# Apply multipole corrections for monopole and dipoles
poissonsolver = MomentCorrectionPoissonSolver(poissonsolver=PoissonSolver(),
                                              moment_corrections=4)

# Ground-state calculation
calc = GPAW(mode='lcao', xc=xc, h=0.3, nbands=216,
            setups={'default':'paw'},
            basis={'default': 'dzp'},
            convergence=convergence,
            poissonsolver=poissonsolver,
            occupations=occupations,
            mixer=mixer,
            maxiter=2000,
            parallel=get_parallel(world),
            symmetry={'point_group': False},
            txt='gs.out')
atoms.set_calculator(calc)
atoms.set_pbc(False)
energy=atoms.get_potential_energy()
#save the ground state energy
calc.write('gs.gpw', mode='all')
