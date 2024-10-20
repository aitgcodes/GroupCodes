import os
import sys
import numpy as np

from gpaw.mpi import world
from gpaw.mpi import SerialCommunicator
from gpaw import GPAW
from ase.io import write
from gpaw import restart
from gpaw.lcaotddft.ksdecomposition import KohnShamDecomposition
from gpaw.tddft.units import as_to_au, au_to_eV
from ase.units import Bohr

sys.path.insert(0, os.path.dirname(__file__))


def generate_cubefiles(pulse, gpw_fpath, ksd_fpath, cube_dpath,
        transitions, basename):
    
    from rhoanalysis.base import BaseCalculator, build_filter

    frequency = pulse.omega0 * au_to_eV
    sigma = pulse.sigma * au_to_eV
    Wlow = frequency - 2 * sigma
    Whigh = frequency + 2 * sigma

    sigma = 0.07

    if transitions == 'all':
        flt1 = None
    elif transitions == 'resonant':
        flt1 = [('w', (Wlow, Whigh))]
    elif transitions == 'nonresonant':
        flt1 = [('w', (0, Wlow)), 'or', ('w', (Whigh, np.inf))]

    # Load ksd
    calc_comm = SerialCommunicator()
    calc = GPAW(gpw_fpath, txt=None, communicator=calc_comm)
    ksd = KohnShamDecomposition(calc, ksd_fpath)
    
    flt_p = build_filter(ksd, flt1)
    eig_n, fermilevel = ksd.get_eig_n(zero_fermilevel=True)
    ia_p = ksd.ia_p[flt_p]
    i_p = ia_p[:, 0]
    a_p = ia_p[:, 1]

    imin = np.min(i_p)
    imax = np.max(i_p)
    amin = np.min(a_p)
    amax = np.max(a_p)

    # load binary file and get calculator
    atoms, calc = restart(gpw_fpath)

    for band in range(imin, amax+1):
        wf = calc.get_pseudo_wave_function(band=band)
        wf_squared = wf**2* Bohr**3

        fname = cube_dpath + '{0}_{1}.cube'.format(basename, band)
        write(fname, atoms, data=wf_squared)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=ExistingPathType)
    parser.add_argument('ksd_fpath', type=ExistingPathType)
    parser.add_argument('pulse_fpath', type=ExistingPathType)
    parser.add_argument('cube_dpath', type=FilePathType)
    parser.add_argument('--transitions', default='all',
                        choices=['all', 'resonant', 'nonresonant'])
    parser.add_argument('basename', type=str)
    args = parser.parse_args()

    from pulse import read_pulse, pulse_time_grid
    pulse = read_pulse(args.pulse_fpath)

    generate_cubefiles(pulse, args.gpw_fpath, args.ksd_fpath, args.cube_dpath,
        args.transitions, args.basename)
