import os
import sys
sys.path.insert(0, os.path.dirname(__file__))


def gs(atoms_fpath,
       out_fpath,
       gpw_fpath,
       nonconv_gpw_fpath,
       mixer=(0.1, 5, 1.0),
       nbands=None,
       fd=0.05):
    from ase.parallel import paropen
    from ase.io import read 
    from gpaw.mpi import world
    from gpaw import GPAW, FermiDirac, Mixer
    from gpaw.poisson import PoissonSolver
    from gpaw.poisson_moment import MomentCorrectionPoissonSolver
    from gpaw import KohnShamConvergenceError
    from gpaw.cluster import Cluster
    from gpaw.utilities.timelimit import TimeLimiter, time_to_seconds
    from parallel_util import get_parallel

    # Check if calculation finished
    if os.path.exists(gpw_fpath):
        return

    # TimeLimiter settings from environment variables
    timelimit = time_to_seconds(os.environ['TIME_LIMIT'])
    timelimit -= time_to_seconds('5:0')
    timestart = float(os.environ['START_TIME'])

    # Calculate
    if not os.path.exists(nonconv_gpw_fpath):
        # Basic settings
        h = 0.3
        v = 6.0

        # Load structure
        #atoms = Cluster(atoms_fpath)
        #atoms.set_pbc(False)
        #atoms.minimal_box(v, h, 32)
        atoms = read(atoms_fpath)
        atoms.set_pbc(False)
        atoms.center(vacuum=6.0)

        # Apply multipole corrections for monopole and dipoles
        poissonsolver = MomentCorrectionPoissonSolver(poissonsolver=PoissonSolver(),
                                              moment_corrections=4)
       
        basis = {'Ag': 'pvalence.dz', 'Au': 'pvalence.dz', 'default': 'dzp'}
        setups = {'Ag': '11', 'default': 'paw'}
 

        # Number of bands
        if nbands is None:
            nbands = {46: 216, 55: 360, 101: 576, 147: 890, 193: 1080,  282: 1290, 337: 1650}.get(len(atoms), None)

        # Calculate
        calc = GPAW(mode='lcao',
                    xc='GLLBSC',
                    basis=basis,
                    convergence={'density': 1e-9},
                    mixer=Mixer(mixer[0], int(mixer[1]), mixer[2]),
                    occupations=FermiDirac(fd),
                    symmetry={'point_group': False},
                    maxiter=2000,
                    h=h,
                    nbands=nbands,
                    setups=setups,
                    poissonsolver=poissonsolver,
                    parallel=get_parallel(world),
                    txt=paropen(out_fpath, 'w')
                    )
        atoms.set_calculator(calc)
    else:  # os.path.exists(nonconv_gpw_fpath):
        # Restart
        calc = GPAW(nonconv_gpw_fpath,
                    parallel=get_parallel(world),
                    txt=paropen(out_fpath, 'a')
                    )
        atoms = calc.get_atoms()

    try:
        tl = TimeLimiter(calc, timelimit=timelimit, timestart=timestart)
        tl.reset('scf')
        atoms.get_potential_energy()
    except KohnShamConvergenceError:
        calc.write(nonconv_gpw_fpath, mode='all')
    else:
        calc.write(gpw_fpath, mode='all')
        if world.rank == 0:
            if os.path.exists(nonconv_gpw_fpath):
                os.remove(nonconv_gpw_fpath)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType, IntOrStrType

    parser = argparse.ArgumentParser()
    parser.add_argument('atoms_fpath', type=ExistingPathType)
    parser.add_argument('out_fpath', type=FilePathType)
    parser.add_argument('gpw_fpath', type=FilePathType)
    parser.add_argument('nonconv_gpw_fpath', type=FilePathType)
    parser.add_argument('--mixer', nargs=3, type=float, default=(0.1, 5, 1.0))
    parser.add_argument('--nbands', type=IntOrStrType)
    parser.add_argument('--fd', type=float, default=0.05)
    args = parser.parse_args()

    gs(args.atoms_fpath,
       args.out_fpath,
       args.gpw_fpath,
       args.nonconv_gpw_fpath,
       nbands=args.nbands,
       mixer=args.mixer,
       fd=args.fd)
