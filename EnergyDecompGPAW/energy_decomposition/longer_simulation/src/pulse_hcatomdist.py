import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def pulse_hcatomdist(pulse, time_t, gpw_fpath, ksd_fpath, pulse_rho_dpath,
                     atomweight_fpath, ofpath, aproj_i):
    from rhoanalysis.hcatomdist import HCAtomDistCalculator

    # Read atomweights
    npz = np.load(atomweight_fpath)
    weight_am = npz['weight_am']
    n_m = npz['n_m']
    d_m = npz['d_m']
    deg_n = npz['deg_n']

    energy_o = np.arange(-5, 1 + 1e-6, 0.01)
    energy_u = np.arange(-1, 5 + 1e-6, 0.01)
    sigma = 0.07

    a = HCAtomDistCalculator(time_t, pulse, gpw_fpath, ksd_fpath,
                             pulse_rho_dpath, weight_am=weight_am,
                             n_m=n_m, d_m=d_m, deg_n=deg_n)
    a.run(energy_o, energy_u, sigma, outfpath=ofpath, aproj_i=aproj_i)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType, IntOrStrType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=ExistingPathType)
    parser.add_argument('ksd_fpath', type=ExistingPathType)
    parser.add_argument('atomweight_fpath', type=ExistingPathType)
    parser.add_argument('pulse_fpath', type=ExistingPathType)
    parser.add_argument('pulse_rho_dpath', type=ExistingPathType)
    parser.add_argument('out_fpath', type=FilePathType)
    parser.add_argument('--time', type=float, nargs='+')
    parser.add_argument('--maxtime', type=float, default=30e3)
    parser.add_argument('--atoms', type=IntOrStrType, nargs='+')
    args = parser.parse_args()

    from pulse import read_pulse, pulse_time_grid
    pulse = read_pulse(args.pulse_fpath)
    time_t = pulse_time_grid(args.time, args.maxtime)

    pulse_hcatomdist(pulse, time_t, args.gpw_fpath, args.ksd_fpath,
                     args.pulse_rho_dpath, args.atomweight_fpath,
                     args.out_fpath, args.atoms)
