import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def pulse_hcdns(pulse, time_t, gpw_fpath, ksd_fpath, pulse_rho_dpath,
                atomweight_fpath, ofpath, elho, energy_limits):
    from rhoanalysis.hcdns import HCDensityCalculator

    # Read state degeneracies
    npz = np.load(atomweight_fpath)
    deg_n = npz['deg_n']

    flt = None
    if energy_limits is not None:
        flt = [(elho[0], energy_limits)]

    a = HCDensityCalculator(time_t, pulse, gpw_fpath, ksd_fpath,
                            pulse_rho_dpath, deg_n=deg_n)
    a.run(elho, outfpath=ofpath, flt=flt)


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
    parser.add_argument('electrons_or_holes', choices=['electrons', 'holes'])
    parser.add_argument('--time', type=float, nargs='+')
    parser.add_argument('--maxtime', type=float, default=30e3)
    parser.add_argument('--energy_limits', type=float, nargs=2, default=None)
    args = parser.parse_args()

    from pulse import read_pulse, pulse_time_grid
    pulse = read_pulse(args.pulse_fpath)
    time_t = pulse_time_grid(args.time, args.maxtime)

    pulse_hcdns(pulse, time_t, args.gpw_fpath, args.ksd_fpath,
                args.pulse_rho_dpath, args.atomweight_fpath,
                args.out_fpath, elho=args.electrons_or_holes,
                energy_limits=args.energy_limits)
