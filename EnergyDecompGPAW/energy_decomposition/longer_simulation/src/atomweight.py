import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def calculate_a_g(gpw_fpath):
    from gpaw.mpi import world
    from gpaw import GPAW
    from gpaw.analyse.wignerseitz import wignerseitz

    calc_comm = world

    calc = GPAW(gpw_fpath, txt=None, communicator=calc_comm,
                parallel={'domain': calc_comm.size})
    calc.initialize_positions()
    atoms = calc.get_atoms()
    gd = calc.density.finegd
    a_g = wignerseitz(gd, atoms)

    big_g = gd.collect(a_g)
    return big_g


def atomweight(gpw_fpath, atomweight_fpath,
               calc_size=8, eig_degeneracy_eps=1e-3,
               eig_limits=[-np.inf, np.inf],
               eig_limits_ref='fermi'):
    import time

    from gpaw.mpi import world
    from gpaw.mpi import SerialCommunicator
    from gpaw import GPAW
    from gpaw.lcaotddft.densitymatrix import get_density
    from parallel_util import get_logger

    # Create communicators
    if calc_size == 1:
        calc_comm = SerialCommunicator()
        loop_comm = world
    elif calc_size == world.size:
        calc_comm = world
        loop_comm = SerialCommunicator()
    else:
        if world.size % calc_size != 0:
            raise RuntimeError('MPI world size is not divisible by '
                               'the calculator communicator size')
        calc_ranks = []
        loop_ranks = []
        for r in range(world.size):
            if r // calc_size == world.rank // calc_size:
                calc_ranks.append(r)
            if r % calc_size == world.rank % calc_size:
                loop_ranks.append(r)
        calc_comm = world.new_communicator(calc_ranks)
        loop_comm = world.new_communicator(loop_ranks)

    log = get_logger(calc_comm, world)
    log('calculate a_g')
    big_a_g = calculate_a_g(gpw_fpath)

    log('start')
    calc = GPAW(gpw_fpath, txt=None, communicator=calc_comm,
                parallel={'domain': calc_comm.size})
    calc.initialize_positions()
    atoms = calc.get_atoms()
    eig_n = calc.get_eigenvalues()
    efermi = calc.get_fermi_level()
    if eig_limits_ref == 'fermi':
        eig_limits_ref = efermi
    eig_limits = np.array(eig_limits)
    eig_limits += eig_limits_ref

    gd = calc.density.finegd
    if calc_comm.rank == 0:
        if world.rank != 0:
            big_a_g = gd.empty(global_array=True, dtype=int)
        loop_comm.broadcast(big_a_g, 0)
    else:
        big_a_g = None
    a_g = gd.zeros(dtype=int)
    gd.distribute(big_a_g, a_g)

    # Find degenerate states
    deg_n = np.zeros(len(eig_n), dtype=int)
    eig_prev = -np.inf
    for n, eig in enumerate(eig_n):
        if eig - eig_prev < eig_degeneracy_eps:
            deg_n[n] = deg_n[n - 1] + 1
        else:
            deg_n[n] = 1
        eig_prev = eig

    Nm = np.sum(deg_n)
    log('Nm', Nm)
    weight_am = np.zeros((len(atoms), Nm))
    n_m = np.zeros(Nm, dtype=int)
    d_m = np.zeros(Nm, dtype=int)
    m = 0
    for n in range(len(eig_n)):
        for d in range(0, deg_n[n]):
            n_m[m] = n
            d_m[m] = d
            m += 1

    u = 0
    wfs = calc.wfs
    C_nM = wfs.kpt_u[u].C_nM

    density_type = 'comp'
    for m in range(loop_comm.rank, Nm, loop_comm.size):
        n = n_m[m]
        eig = eig_n[n]
        if eig < eig_limits[0] or eig >= eig_limits[-1]:
            weight_am[:, m] = np.nan
            continue

        d = d_m[m]
        log('n', n, 'd', d)
        t0 = time.time()
        rho_MM = np.outer(C_nM[n], C_nM[n - d])
        rho_MM = 0.5 * (rho_MM + rho_MM.T)
        rho_g = get_density(rho_MM, calc.wfs, calc.density,
                            density_type, 0)
        for a, atom in enumerate(atoms):
            weight = gd.integrate(np.where(a_g == a, rho_g, 0.0))
            weight_am[a, m] = weight
        if d == 0:
            log('sum', np.sum(weight_am[:, m], axis=0))
        t1 = time.time()
        log('time', '%.1f s' % (t1 - t0))

    loop_comm.sum(weight_am, 0)
    if world.rank == 0:
        np.savez_compressed(atomweight_fpath,
                            weight_am=weight_am,
                            n_m=n_m, d_m=d_m,
                            deg_n=deg_n)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType, FloatOrStrType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=ExistingPathType)
    parser.add_argument('atomweight_fpath', type=FilePathType)
    parser.add_argument('--calc_size', type=int, default=8)
    parser.add_argument('--eig_limits', type=float, nargs=2,
                        default=[-np.inf, np.inf])
    parser.add_argument('--eig_limits_ref', type=FloatOrStrType,
                        default='fermi')
    args = parser.parse_args()

    atomweight(args.gpw_fpath, args.atomweight_fpath,
               calc_size=args.calc_size,
               eig_limits=args.eig_limits,
               eig_limits_ref=args.eig_limits_ref)
