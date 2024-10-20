import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def generate_cubefiles(gpw_fpath, out_dpath, basename, nbands_min, nbands_max):
    from gpaw.mpi import world
    from gpaw import GPAW
    from ase.io import write
    from gpaw import restart
    import numpy as np
    from ase.units import Bohr
    from gpaw.utilities.ps2ae import PS2AE

    # load binary file and get calculator
    atoms, calc = restart(gpw_fpath)

    # loop over all wfs and write their cube files
    #nbands = calc.get_number_of_bands()

    for band in range(nbands_min, nbands_max):
        wf = calc.get_pseudo_wave_function(band=band)
        wf_squared = wf**2* Bohr**3

        fname = out_dpath + '{0}_{1}.cube'.format(basename, band)
        write(fname, atoms, data=wf_squared)

        #gd = calc.wfs.gd

        #X = gd.coords(0) * Bohr
        #Y = gd.coords(1) * Bohr
        #Z = gd.coords(2) * Bohr


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=ExistingPathType)
    parser.add_argument('out_dpath', type=FilePathType)
    parser.add_argument('basename', type=str)
    parser.add_argument('--nbands_min', type=int, default=0)
    parser.add_argument('--nbands_max', type=int, default=0)
    args = parser.parse_args()

    generate_cubefiles(args.gpw_fpath, args.out_dpath, args.basename, args.nbands_min, args.nbands_max)
