import os
import sys
sys.path.insert(0, os.path.dirname(__file__))


def rlx(atoms_fpath,
        rlx_traj_fpath,
        rlx_old_traj_fpath,
        rlx_new_traj_fpath,
        rlx_out_fpath,
        rlx_log_fpath,
        nbands=None, h=0.2, v=6.0):
    from ase.io import read
    from ase.parallel import paropen
    from ase.optimize import BFGS
    from gpaw.mpi import world
    from gpaw import GPAW, FermiDirac, PoissonSolver, Mixer
    from gpaw.cluster import Cluster
    from parallel_util import get_parallel

    # Check if we can restart
    do_restart = (os.path.exists(rlx_traj_fpath) or
                  os.path.exists(rlx_old_traj_fpath))

    # Load structure
    if do_restart:
        if world.rank == 0:
            if not os.path.exists(rlx_traj_fpath):
                os.system(f'ase convert '
                          f'{rlx_old_traj_fpath} '
                          f'{rlx_new_traj_fpath} '
                          f'{rlx_traj_fpath}')
        world.barrier()
        atoms = read(rlx_traj_fpath)
    else:
        atoms = Cluster(atoms_fpath)
        atoms.minimal_box(v, h, 32)
        atoms.set_pbc(False)

    # PoissonSolver
    poissonsolver = PoissonSolver(eps=1e-16, remove_moment=4)

    # Setups
    setups = 'my'

    # Calculate
    calc = GPAW(mode='lcao',
                xc='PBE',
                basis='PBE.dzp',
                convergence={'density': 1e-6},
                mixer=Mixer(0.1, 5, 50),
                occupations=FermiDirac(0.05),
                symmetry={'point_group': False},
                h=h,
                nbands=nbands,
                setups=setups,
                poissonsolver=poissonsolver,
                parallel=get_parallel(world),
                txt=paropen(rlx_out_fpath, 'a' if do_restart else 'w')
                )
    atoms.set_calculator(calc)

    if do_restart:
        dyn = BFGS(atoms, logfile=rlx_log_fpath,
                   trajectory=rlx_new_traj_fpath)
        dyn.replay_trajectory(rlx_traj_fpath)
        if world.rank == 0:
            os.rename(rlx_traj_fpath, rlx_old_traj_fpath)
        world.barrier()
    else:
        dyn = BFGS(atoms, logfile=rlx_log_fpath,
                   trajectory=rlx_traj_fpath)
    dyn.run(fmax=0.01)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType, IntOrStrType

    parser = argparse.ArgumentParser()
    parser.add_argument('atoms_fpath', type=ExistingPathType)
    parser.add_argument('rlx_traj_fpath', type=FilePathType)
    parser.add_argument('rlx_old_traj_fpath', type=FilePathType)
    parser.add_argument('rlx_new_traj_fpath', type=FilePathType)
    parser.add_argument('rlx_out_fpath', type=FilePathType)
    parser.add_argument('rlx_log_fpath', type=FilePathType)
    parser.add_argument('--nbands', type=IntOrStrType, default='103%')
    args = parser.parse_args()

    rlx(args.atoms_fpath,
        args.rlx_traj_fpath,
        args.rlx_old_traj_fpath,
        args.rlx_new_traj_fpath,
        args.rlx_out_fpath,
        args.rlx_log_fpath,
        nbands=args.nbands)
