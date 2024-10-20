import os
import sys
import numpy as np

from gpaw.mpi import world
from gpaw.mpi import SerialCommunicator
from gpaw import GPAW
from gpaw.lcaotddft.ksdecomposition import KohnShamDecomposition
from gpaw.tddft.units import as_to_au, au_to_eV

sys.path.insert(0, os.path.dirname(__file__))


def weighted_energy(pulse, time_t, gpw_fpath, ksd_fpath, pulse_rho_dpath,
        ofpath, transitions, wi_full_fpath, wa_full_fpath, wi_np_fpath, wa_np_fpath):
    
    from rhoanalysis.base import BaseCalculator, build_filter
    from rhoanalysis.energy_Eia import EnergyCalculator
    from rhoanalysis.wtrans_energy import WTECalculator
    from read_weights_wfs import ks_weights

    frequency = pulse.omega0 * au_to_eV
    sigma = pulse.sigma * au_to_eV
    Wlow = frequency - 2 * sigma
    Whigh = frequency + 2 * sigma

    energy_o = np.arange(-5, 1 + 1e-6, 0.01)
    energy_u = np.arange(-1, 5 + 1e-6, 0.01)
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

    weights_wfs=ks_weights(wi_full_fpath,wa_full_fpath,wi_np_fpath,wa_np_fpath,i_p,a_p)

    # The flt has to be rebuilt as it gets empty in the previous step
    if transitions == 'all':
        flt2 = None
    elif transitions == 'resonant':
        flt2 = [('w', (Wlow, Whigh))]
    elif transitions == 'nonresonant':
        flt2 = [('w', (0, Wlow)), 'or', ('w', (Whigh, np.inf))]


    a = EnergyCalculator(time_t, pulse, gpw_fpath, ksd_fpath, pulse_rho_dpath)
    E_ia = a.run(flt=flt2)

    # The flt has to be rebuilt as it gets empty in the previous step
    if transitions == 'all':
        flt3 = None
    elif transitions == 'resonant':
        flt3 = [('w', (Wlow, Whigh))]
    elif transitions == 'nonresonant':
        flt3 = [('w', (0, Wlow)), 'or', ('w', (Whigh, np.inf))]


    b = WTECalculator(time_t,pulse,gpw_fpath,ksd_fpath,pulse_rho_dpath)
    b.run(weights_wfs, E_ia, outfpath=ofpath, flt=flt3)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=ExistingPathType)
    parser.add_argument('ksd_fpath', type=ExistingPathType)
    parser.add_argument('pulse_fpath', type=ExistingPathType)
    parser.add_argument('pulse_rho_dpath', type=ExistingPathType)
    parser.add_argument('wi_full_fpath', type=FilePathType)
    parser.add_argument('wa_full_fpath', type=FilePathType)
    parser.add_argument('wi_np_fpath', type=FilePathType)
    parser.add_argument('wa_np_fpath', type=FilePathType)
    parser.add_argument('out_fpath', type=FilePathType)
    parser.add_argument('--time', type=float, nargs='+')
    parser.add_argument('--maxtime', type=float, default=30e3)
    parser.add_argument('--transitions', default='all',
                        choices=['all', 'resonant', 'nonresonant'])
    args = parser.parse_args()

    from pulse import read_pulse, pulse_time_grid
    pulse = read_pulse(args.pulse_fpath)
    time_t = pulse_time_grid(args.time, args.maxtime)

    weighted_energy(pulse,time_t,args.gpw_fpath,args.ksd_fpath,args.pulse_rho_dpath,
                 args.out_fpath,args.transitions,args.wi_full_fpath,args.wa_full_fpath,args.wi_np_fpath,args.wa_np_fpath)
