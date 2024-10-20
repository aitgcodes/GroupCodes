import os
import sys
import numpy as np
import subprocess

from gpaw.mpi import world
from gpaw.mpi import SerialCommunicator
from gpaw import GPAW
from ase.io import write
from gpaw import restart
from gpaw.lcaotddft.ksdecomposition import KohnShamDecomposition
from gpaw.tddft.units import as_to_au, au_to_eV
from ase.units import Bohr
sys.path.insert(0, os.path.dirname(__file__))

def ks_weights(pulse, gpw_fpath, ksd_fpath, cube_dpath, 
            basename, transitions, occ_weigth_out_fpath,
            unocc_weigth_out_fpath, lx, ly, lz):

    from rhoanalysis.base import BaseCalculator, build_filter

    frequency = pulse.omega0 * au_to_eV
    sigma = pulse.sigma * au_to_eV
    Wlow = frequency - 2 * sigma
    Whigh = frequency + 2 * sigma

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


    f1 = open(occ_weigth_out_fpath, "w")
    f2 = open(unocc_weigth_out_fpath, "w")

    for i in range(imin, imax+1):
        out = subprocess.Popen(["python3","src/cubeintegrator.py",cube_dpath+basename+"_"+str(i)+".cube", lx,ly,lz],stdout=subprocess.PIPE,universal_newlines=True)
        out.wait()  # wait until the script has finished
        stdout_data, stderr_data = out.communicate()
        f1.writelines(["State  "+str(i)+"  ",stdout_data])

    for a in range(amin, amax+1):
        out = subprocess.Popen(["python3","src/cubeintegrator.py",cube_dpath+basename+"_"+str(a)+".cube", lx,ly,lz],stdout=subprocess.PIPE,universal_newlines=True)
        out.wait()  # wait until the script has finished
        stdout_data, stderr_data = out.communicate()
        f2.writelines(["State  "+str(a)+"  ",stdout_data])

    f1.close()
    f2.close()

if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=ExistingPathType)
    parser.add_argument('ksd_fpath', type=ExistingPathType)
    parser.add_argument('pulse_fpath', type=ExistingPathType)
    parser.add_argument('cube_dpath', type=ExistingPathType)
    parser.add_argument('basename', type=str)
    parser.add_argument('--transitions', default='all',
                        choices=['all', 'resonant', 'nonresonant'])
    parser.add_argument('occ_weigth_out_fpath', type=FilePathType)
    parser.add_argument('unocc_weigth_out_fpath', type=FilePathType)
    parser.add_argument('--lx', type=str, default='-')
    parser.add_argument('--ly', type=str, default='-')
    parser.add_argument('--lz', type=str, default='-')

    args = parser.parse_args()

    from pulse import read_pulse, pulse_time_grid
    pulse = read_pulse(args.pulse_fpath)

    ks_weights(pulse, args.gpw_fpath, args.ksd_fpath, args.cube_dpath, 
            args.basename, args.transitions,args.occ_weigth_out_fpath,
            args.unocc_weigth_out_fpath, args.lx, args.ly, args.lz)
