import os
import numpy as np
from ase.io import write
from gpaw import GPAW
from gpaw.tddft.units import au_to_eV
from gpaw.lcaotddft import LCAOTDDFT
from gpaw.lcaotddft.densitymatrix import DensityMatrix
from gpaw.lcaotddft.frequencydensitymatrix import FrequencyDensityMatrix
from gpaw.lcaotddft.ksdecomposition import KohnShamDecomposition
#from ksdecomposition import KohnShamDecomposition
from pathlib import Path


def read_rho(rhopath, time, tag):
    fpath = os.path.join(rhopath, 't%09.1f%s.npy' % (time, tag))
    if not os.path.exists(fpath):
            raise RuntimeError('File missing: %s' % fpath)
    rho_p = np.load(fpath)
    return rho_p

def ind_den(gpw_fpath, ksd_fpath, pulse_rho_dpath, ofpath, time):
    calc= GPAW(gpw_fpath, txt=None)
    calc.initialize_positions()  # Initialize in order to calculate density
    dmat = DensityMatrix(calc)
    ksd = KohnShamDecomposition(calc, ksd_fpath)

    Path(ofpath).mkdir(parents=True, exist_ok=True)

    rho_MM_0 = read_rho(pulse_rho_dpath, 0, "")
    rho_MM_t = read_rho(pulse_rho_dpath, time, "")
    rho_up = rho_MM_t - rho_MM_0
    rho_g = ksd.get_density(calc.wfs, rho_up.imag)
    t_fs = time/1000 # time from as to fs
    write(f'{ofpath}/ind_{t_fs:.2f}fs.cube', calc.atoms, data=rho_g)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=ExistingPathType)
    parser.add_argument('ksd_fpath', type=ExistingPathType)
    parser.add_argument('pulse_rho_dpath', type=ExistingPathType)
    parser.add_argument('out_fpath', type=FilePathType)
    parser.add_argument('--time', type=float)

    args = parser.parse_args()

    ind_den(args.gpw_fpath, args.ksd_fpath, args.pulse_rho_dpath,
            args.out_fpath, args.time)
