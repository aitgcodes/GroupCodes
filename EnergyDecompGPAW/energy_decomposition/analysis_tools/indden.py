import os
import sys
sys.path.insert(0, os.path.dirname(__file__))

from pathlib import Path
import numpy as np


def ind_den(gpw_fpath, ksd_fpath, rho_fpath, out_fpath, time):

    from gpaw import GPAW
    from gpaw.tddft.units import au_to_eV
    from gpaw.lcaotddft.ksdecomposition import KohnShamDecomposition
    from gpaw.lcaotddft.densitymatrix import DensityMatrix
    from gpaw.lcaotddft.frequencydensitymatrix import FrequencyDensityMatrix

    # Load the objects
    calc = GPAW(gpw_fpath, txt=None)
    calc.initialize_positions()  # Initialize in order to calculate density
    dmat = DensityMatrix(calc)
    ksd = KohnShamDecomposition(calc, ksd_fpath)

    # Load the density matrix at a given time
    rho_p = np.load(rho_fpath)
    rho_g = dmat.get_density([rho_p.imag])
    write(out_fpath + f'ind_{time:.2f}.cube', calc.atoms, data=rho_g)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType, IntOrStrType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=FilePathType)
    parser.add_argument('ksd_fpath', type=FilePathType)
    parser.add_argument('rho_fpath', type=FilePathType)
    parser.add_argument('out_fpath', type=FilePathType)
    parser.add_argument('--time', type=float)

    args = parser.parse_args()

    ind_den(args.gpw_fpath,
            args.ksd_fpath, 
            args.rho_fpath, 
            args.out_fpath, 
            args.time)
