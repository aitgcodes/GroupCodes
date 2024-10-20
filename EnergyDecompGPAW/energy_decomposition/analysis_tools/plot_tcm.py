import os
import sys
sys.path.insert(0, os.path.dirname(__file__))

from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt


def plot_tcm(gpw_fpath, ksd_fpath, tcm_fpath, out_fpath, freq):

    from gpaw import GPAW
    from gpaw.tddft.units import au_to_eV
    from gpaw.lcaotddft.ksdecomposition import KohnShamDecomposition
    from gpaw.lcaotddft.densitymatrix import DensityMatrix
    from gpaw.lcaotddft.frequencydensitymatrix import FrequencyDensityMatrix
    from gpaw.lcaotddft.tcm import TCMPlotter

    from tcm import TCMPlotter

    # Load the objects
    calc = GPAW(gpw_fpath, txt=None)
    ksd = KohnShamDecomposition(calc, ksd_fpath)

    # Load the time-dependent tcm
    tcm_en = np.load(tcm_fpath)
    time = tcm_en["time"]/1000     # Time in fs
    energy_o = tcm_en["energy_o"]
    energy_u = tcm_en["energy_u"]
    sigma = tcm_en["sigma"]
    tcm_ou = tcm_en["tcm_ou"]

    #----------------------------------------

    wlow = freq - 8*sigma   # lower frequency limit, below which nonresonant plasmonic transitions occur
    whigh = freq + 8*sigma  # upper frequency limit, above which nonresonant plasmonic transitions occur

    plt.figure(figsize=(8, 8))
    plt.clf()

    plotter = TCMPlotter(ksd, energy_o, energy_u, sigma)
    plotter.plot_TCM(tcm_ou)
    plotter.plot_TCM_diagonal(freq, color='k')
    plotter.plot_TCM_diagonal(wlow, color='k', linestyle=':')
    plotter.plot_TCM_diagonal(whigh, color='k', linestyle=':')

    # Save the plot
    plt.savefig(out_fpath + f'tcm_{freq:.2f}_{time:.2f}.png')

if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType, IntOrStrType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=FilePathType)
    parser.add_argument('ksd_fpath', type=FilePathType)
    parser.add_argument('tcm_fpath', type=FilePathType)
    parser.add_argument('out_fpath', type=FilePathType)
    parser.add_argument('--freq', type=float)

    args = parser.parse_args()

    plot_tcm(args.gpw_fpath,
            args.ksd_fpath, 
            args.tcm_fpath, 
            args.out_fpath, 
            args.freq)
