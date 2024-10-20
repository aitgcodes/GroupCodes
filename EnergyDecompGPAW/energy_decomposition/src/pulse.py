import os
import sys
import pickle
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def create_pulse(frequency):
    from gpaw.lcaotddft.laser import GaussianPulse
    from gpaw.tddft.units import fs_to_au, au_to_eV

    # Pulse
    fwhm = 5.0  # fs
    tau = fwhm / (2 * np.sqrt(2 * np.log(2)))
    sigma = 1 / (tau * fs_to_au) * au_to_eV  # eV
    strength = 1e-5
    t0 = 10e3
    sincos = 'sin'
    return GaussianPulse(strength, t0, frequency, sigma, sincos)


def convolution_freq_grid(pulse, freq_step):
    from gpaw.tddft.units import au_to_eV

    frequency = pulse.omega0 * au_to_eV
    sigma = pulse.sigma * au_to_eV
    buf = 4 * sigma

    # Create frequency grid covering multiplets of `freq_step` in
    # [freq - buf, freq + buf] (flooring and ceiling)
    freq_w = np.arange(0, frequency + 2 * buf, freq_step)
    flt_w = np.logical_and(frequency - buf < freq_w, freq_w < frequency + buf)
    flt_w[:-1] += flt_w[1:]
    flt_w[1:] += flt_w[:-1]
    freq_w = freq_w[flt_w]
    return freq_w


def pulse_time_grid(time, maxtime):
    if time is None:
        dtime = 20.0
        time_t = np.arange(0, maxtime + 0.5 * dtime, dtime)
    else:
        time_t = np.array(time)
    return time_t


def write_pulse(pulse, fpath):
    with open(fpath, 'wb') as f:
        pickle.dump(pulse.todict(), f, pickle.HIGHEST_PROTOCOL)


def read_pulse(fpath):
    from gpaw.lcaotddft.laser import create_laser

    with open(fpath, 'rb') as f:
        dct = pickle.load(f)
    return create_laser(**dct)


if __name__ == '__main__':
    import argparse
    from argparse_util import FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('pulse_fpath', type=FilePathType)
    parser.add_argument('frequency', type=float)
    args = parser.parse_args()

    pulse = create_pulse(args.frequency)
    write_pulse(pulse, args.pulse_fpath)
