import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def pulse_dm(pulse, freq_w, time_t, kick_dm_fpath, pulse_dm_fpath):
    from gpaw.tddft.spectrum import read_dipole_moment_file
    from gpaw.tddft.units import eV_to_au, as_to_au
    from convolution import TimeConvolver

    # Read dipole moment file
    kick_j, time_T, norm_T, dm_Tv = read_dipole_moment_file(kick_dm_fpath)
    if len(kick_j) != 1:
        raise RuntimeError('there must be exactly one kick')
    kick_v = kick_j[0]['strength_v']
    dm_Tv /= np.sqrt(np.sum(kick_v**2))

    # Do convolution in time domain
    c = TimeConvolver(pulse, time_T, dm_Tv)
    pulsedm_tv = c.convolve(time_t, units='as')

    # Do convolution in frequency domain
    pulsedmft_tv = c.fourier(freq_w * eV_to_au).convolve(time_t, units='as')

    # Save data
    autime_t = time_t * as_to_au
    pulse_t = pulse.strength(autime_t)
    np.savetxt(pulse_dm_fpath, np.concatenate((autime_t[:, np.newaxis],
                                               pulse_t[:, np.newaxis],
                                               pulsedmft_tv,
                                               pulsedm_tv), axis=1))


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('kick_dm_fpath', type=ExistingPathType)
    parser.add_argument('pulse_fpath', type=ExistingPathType)
    parser.add_argument('pulse_dm_fpath', type=FilePathType)
    parser.add_argument('--freq_step', type=float, default=0.05)
    parser.add_argument('--time', type=float, nargs='+')
    parser.add_argument('--maxtime', type=float, default=30e3)
    args = parser.parse_args()

    from pulse import read_pulse, convolution_freq_grid, pulse_time_grid
    pulse = read_pulse(args.pulse_fpath)
    freq_w = convolution_freq_grid(pulse, args.freq_step)
    time_t = pulse_time_grid(args.time, args.maxtime)

    pulse_dm(pulse, freq_w, time_t, args.kick_dm_fpath, args.pulse_dm_fpath)
