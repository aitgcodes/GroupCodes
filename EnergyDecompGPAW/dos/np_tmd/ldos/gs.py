import os
import sys
from ase.io import read
from ase.parallel import paropen
from ase.io import read
from gpaw.mpi import world
from gpaw import GPAW, FermiDirac, PoissonSolver, Mixer
from gpaw import KohnShamConvergenceError
from gpaw.cluster import Cluster
from gpaw.utilities.timelimit import TimeLimiter, time_to_seconds
from parallel_util import get_parallel
from gpaw.poisson_moment import MomentCorrectionPoissonSolver

# Read the structure from the xyz file
atoms = read('rlx-ico-Ag55-MoSe2_4-geom1.xyz')
atoms.set_pbc(False)
atoms.center(vacuum=6.0)

# Apply multipole corrections for monopole and dipoles
poissonsolver = MomentCorrectionPoissonSolver(poissonsolver=PoissonSolver(),
                                              moment_corrections=4)

mixer=(0.1, 5, 1.0)
fd=0.05
nbands = {46: 216, 55: 360, 101: 576, 147: 890, 193: 1080,  282: 1290, 337: 1650}.get(len(atoms), None)
h=0.3

basis = {'Ag': 'pvalence.dz', 'Au': 'pvalence.dz', 'default': 'dzp'}
setups={'Ag': '11', 'default': 'paw'}

# Calculate
calc = GPAW(mode='lcao',
            xc='GLLBSC',
            basis=basis,
            convergence={'density': 1e-9},
            mixer=Mixer(mixer[0], int(mixer[1]), mixer[2]),
            occupations=FermiDirac(fd),
            symmetry={'point_group': False},
            maxiter=1000,
            h=h,
            nbands=nbands,
            setups=setups,
            poissonsolver=poissonsolver,
            parallel=get_parallel(world),
            txt='gs.out')

atoms.set_calculator(calc)
energy=atoms.get_potential_energy()
#save the ground state energy
calc.write('gs.gpw', mode='all')
