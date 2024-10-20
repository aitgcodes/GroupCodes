from ase.io import read
from gpaw import GPAW, FermiDirac, Mixer, PoissonSolver
from ase.optimize import BFGS
from gpaw.mpi import world
from parallel_util import get_parallel


# Read the structure from the xyz file
atoms = read('unrlx-mose2_4.xyz')
atoms.center(vacuum=6.0)
atoms.set_pbc(False)

# PoissonSolver
poissonsolver = PoissonSolver(remove_moment=4)
h=0.2
setups={'default': 'paw'}
nbands=220

# Calculate
calc = GPAW(mode='lcao',
            xc='PBE',
            basis='dzp',
            convergence={'density': 1e-6},
            mixer=Mixer(0.1, 5, 1),
            occupations=FermiDirac(0.05),
            symmetry={'point_group': False},
            h=h,
            nbands=nbands,
            setups=setups,
            poissonsolver=poissonsolver,
            parallel=get_parallel(world),
            txt='opt.out')
                
atoms.set_calculator(calc)

dyn = BFGS(atoms,trajectory='opt.traj', logfile='opt.log')
dyn.run(fmax=0.01)
