import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def pulse_energy(pulse, time_t, gpw_fpath, ksd_fpath, pulse_rho_dpath,
                 ofpath, transitions):
    from gpaw.tddft.units import au_to_eV
    from rhoanalysis.energy import EnergyCalculator

    frequency = pulse.omega0 * au_to_eV
    sigma = pulse.sigma * au_to_eV
    Wlow = frequency - 2 * sigma
    Whigh = frequency + 2 * sigma

    if transitions == 'all':
        flt = None
    elif transitions == 'resonant':
        flt = [('w', (Wlow, Whigh))]
    elif transitions == 'nonresonant':
        flt = [('w', (0, Wlow)), 'or', ('w', (Whigh, np.inf))]

    a = EnergyCalculator(time_t, pulse, gpw_fpath, ksd_fpath, pulse_rho_dpath)
    a.run(outfpath=ofpath, flt=flt)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=ExistingPathType)
    parser.add_argument('ksd_fpath', type=ExistingPathType)
    parser.add_argument('pulse_fpath', type=ExistingPathType)
    parser.add_argument('pulse_rho_dpath', type=ExistingPathType)
    parser.add_argument('out_fpath', type=FilePathType)
    parser.add_argument('--time', type=float, nargs='+')
    parser.add_argument('--maxtime', type=float, default=30e3)
    parser.add_argument('--transitions', default='all',
                        choices=['all', 'resonant', 'nonresonant'])
    args = parser.parse_args()

    from pulse import read_pulse, pulse_time_grid
    pulse = read_pulse(args.pulse_fpath)
    time_t = pulse_time_grid(args.time, args.maxtime)

    pulse_energy(pulse, time_t, args.gpw_fpath, args.ksd_fpath,
                 args.pulse_rho_dpath, args.out_fpath,
                 transitions=args.transitions)
