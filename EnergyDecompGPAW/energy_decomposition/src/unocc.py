import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def unocc(gs_gpw_fpath,
          unocc_out_fpath,
          unocc_gpw_fpath,
          eig_fpath):
    from ase.parallel import paropen
    from gpaw.mpi import world
    from gpaw import GPAW
    from parallel_util import get_parallel

    # Check if calculation finished
    if os.path.exists(unocc_gpw_fpath):
        return

    # Restart
    calc = GPAW(gs_gpw_fpath,
                nbands='nao',
                fixdensity=True,
                convergence={'bands': 'all'},
                parallel=get_parallel(world),
                txt=paropen(unocc_out_fpath, 'w')
                )
    atoms = calc.get_atoms()
    atoms.get_potential_energy()
    calc.write(unocc_gpw_fpath, mode='all')

    # Write eigenvalues
    fermie = calc.get_fermi_level()
    eig_n = calc.get_eigenvalues()
    occ_n = calc.get_occupation_numbers()
    with paropen(eig_fpath, 'w') as f:
        f.write('# Fermi level: %.8f\n' % fermie)
        f.write('# %14s %16s\n' % ('Eigenvalue', 'Occupation'))
        np.savetxt(f, np.stack((eig_n, occ_n)).T, fmt='%16.8f')


def ksd(unocc_gpw_fpath, ksd_out_fpath, ksd_fpath):
    from gpaw import GPAW
    from gpaw.lcaotddft.ksdecomposition import KohnShamDecomposition

    calc = GPAW(unocc_gpw_fpath, txt=ksd_out_fpath)

    # Construct KS electron-hole basis
    ksd = KohnShamDecomposition(calc)
    ksd.initialize(calc)
    ksd.write(ksd_fpath)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('gs_gpw_fpath', type=ExistingPathType)
    parser.add_argument('unocc_out_fpath', type=FilePathType)
    parser.add_argument('unocc_gpw_fpath', type=FilePathType)
    parser.add_argument('eig_fpath', type=FilePathType)
    parser.add_argument('ksd_out_fpath', type=FilePathType)
    parser.add_argument('ksd_fpath', type=FilePathType)
    args = parser.parse_args()

    unocc(args.gs_gpw_fpath,
          args.unocc_out_fpath,
          args.unocc_gpw_fpath,
          args.eig_fpath)
    ksd(args.unocc_gpw_fpath,
        args.ksd_out_fpath,
        args.ksd_fpath)
