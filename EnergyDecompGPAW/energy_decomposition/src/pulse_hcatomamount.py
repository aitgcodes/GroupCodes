import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def pulse_hcatomamount(pulse, time_t, gpw_fpath, ksd_fpath, pulse_rho_dpath,
                       atomweight_fpath, ofpath, holes, electrons):
    from rhoanalysis.hcatomamount import HCAtomAmountCalculator

    # Read atomweights
    npz = np.load(atomweight_fpath)
    weight_am = npz['weight_am']
    n_m = npz['n_m']
    d_m = npz['d_m']
    deg_n = npz['deg_n']

    flt = None
    if holes is not None:
        flt = [('i', holes)]
        if electrons is not None:
            flt += ['and', ('a', electrons)]
    else:
        if electrons is not None:
            flt = [('a', electrons)]

    a = HCAtomAmountCalculator(time_t, pulse, gpw_fpath, ksd_fpath,
                               pulse_rho_dpath, weight_am=weight_am,
                               n_m=n_m, d_m=d_m, deg_n=deg_n)
    a.run(outfpath=ofpath, flt=flt)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=ExistingPathType)
    parser.add_argument('ksd_fpath', type=ExistingPathType)
    parser.add_argument('atomweight_fpath', type=ExistingPathType)
    parser.add_argument('pulse_fpath', type=ExistingPathType)
    parser.add_argument('pulse_rho_dpath', type=ExistingPathType)
    parser.add_argument('out_fpath', type=FilePathType)
    parser.add_argument('--time', type=float, nargs='+')
    parser.add_argument('--maxtime', type=float, default=30e3)
    parser.add_argument('--holes', type=float, nargs=2, default=None)
    parser.add_argument('--electrons', type=float, nargs=2, default=None)
    args = parser.parse_args()

    from pulse import read_pulse, pulse_time_grid
    pulse = read_pulse(args.pulse_fpath)
    time_t = pulse_time_grid(args.time, args.maxtime)

    pulse_hcatomamount(pulse, time_t, args.gpw_fpath, args.ksd_fpath,
                       args.pulse_rho_dpath, args.atomweight_fpath,
                       args.out_fpath, args.holes, args.electrons)
