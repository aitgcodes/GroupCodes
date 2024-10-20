import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def dos(gpw_fpath, dos_fpath, energy_e):
    from gpaw.mpi import world
    from gpaw import GPAW
    from gpaw.lcaotddft.ksdecomposition import gauss_ij

    calc = GPAW(gpw_fpath, txt=None)
    eig_n = calc.get_eigenvalues()
    efermi = calc.get_fermi_level()

    # DOS
    width = 0.07
    G_en = gauss_ij(energy_e, eig_n - efermi, width)
    dos_e = np.sum(G_en, axis=1)

    if world.rank == 0:
        with open(dos_fpath, 'w') as f:
            f.write(f'# Density of states wrt Fermi level')
            f.write(f'# Gaussian folding, width = {width:.4f} eV\n')
            np.savetxt(f, np.stack((energy_e, dos_e)).T)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=ExistingPathType)
    parser.add_argument('dos_fpath', type=FilePathType)
    args = parser.parse_args()

    # Energy grid
    emin = -20.0
    emax = 20.0
    de = 0.01
    energy_e = np.arange(emin, emax + 0.5 * de, de)

    dos(args.gpw_fpath, args.dos_fpath, energy_e)
