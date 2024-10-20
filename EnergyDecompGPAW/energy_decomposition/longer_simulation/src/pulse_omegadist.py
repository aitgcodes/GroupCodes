import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def pulse_omegadist(pulse, time_t, gpw_fpath, ksd_fpath, pulse_rho_dpath,
                    ofpath, weight):
    from rhoanalysis.omegadist import OmegaDistCalculator

    energy_e = np.arange(0, 9 + 1e-6, 0.01)
    sigma = 0.07

    a = OmegaDistCalculator(time_t, pulse, gpw_fpath, ksd_fpath,
                            pulse_rho_dpath)
    a.run(energy_e, sigma, outfpath=ofpath, weight=weight)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=ExistingPathType)
    parser.add_argument('ksd_fpath', type=ExistingPathType)
    parser.add_argument('pulse_fpath', type=ExistingPathType)
    parser.add_argument('pulse_rho_dpath', type=ExistingPathType)
    parser.add_argument('out_fpath', type=FilePathType)
    parser.add_argument('quantity', choices=['occupation',
                                             'energy',
                                             'coulombenergy',
                                             'rate_absorption',
                                             'rate_interaction'])
    parser.add_argument('--time', type=float, nargs='+')
    parser.add_argument('--maxtime', type=float, default=30e3)
    args = parser.parse_args()

    from pulse import read_pulse, pulse_time_grid
    pulse = read_pulse(args.pulse_fpath)
    time_t = pulse_time_grid(args.time, args.maxtime)

    pulse_omegadist(pulse, time_t, args.gpw_fpath, args.ksd_fpath,
                    args.pulse_rho_dpath, args.out_fpath,
                    weight=args.quantity)
