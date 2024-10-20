import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def pulse_tcm(pulse, time_t, gpw_fpath, ksd_fpath, pulse_rho_dpath,
              odpath, weight):
    from rhoanalysis.tcm import TCMCalculator

    energy_o = np.arange(-5, 0.5 + 1e-6, 0.01)
    energy_u = np.arange(-0.5, 5 + 1e-6, 0.01)
    sigma = 0.07

    a = TCMCalculator(time_t, pulse, gpw_fpath, ksd_fpath, pulse_rho_dpath)
    a.run(energy_o, energy_u, sigma, outdpath=odpath, weight=weight)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, DirectoryPathType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=ExistingPathType)
    parser.add_argument('ksd_fpath', type=ExistingPathType)
    parser.add_argument('pulse_fpath', type=ExistingPathType)
    parser.add_argument('pulse_rho_dpath', type=ExistingPathType)
    parser.add_argument('out_dpath', type=DirectoryPathType)
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

    pulse_tcm(pulse, time_t, args.gpw_fpath, args.ksd_fpath,
              args.pulse_rho_dpath, args.out_dpath,
              weight=args.quantity)
