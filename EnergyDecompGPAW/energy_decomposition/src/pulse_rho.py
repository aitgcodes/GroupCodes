import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def pulse_rho(pulse, freq_w, time_t, gpw_fpath, ksd_fpath, frho_dpath,
              pulse_rho_dpath, kickstr=1e-5):
    from gpaw.mpi import world
    from gpaw.mpi import SerialCommunicator
    from gpaw import GPAW
    from gpaw.lcaotddft.ksdecomposition import KohnShamDecomposition
    from gpaw.tddft.units import eV_to_au
    from parallel_util import get_logger
    from convolution import FourierConvolver

    calc_comm = SerialCommunicator()
    ksd_comm = world
    log = get_logger(ksd_comm, ksd_comm)

    # Load ksd
    calc = GPAW(gpw_fpath, txt=None, communicator=calc_comm)
    ksd = KohnShamDecomposition(calc, ksd_fpath)
    ksd.distribute(ksd_comm)

    # Read density matrix
    reim_r = ['Re', 'Im']
    rho_wrq = np.zeros((len(freq_w), len(reim_r), len(ksd.w_q)),
                       dtype=complex)
    for w, freq in enumerate(freq_w):
        for r, reim in enumerate(reim_r):
            fpath = os.path.join(frho_dpath, 'w%05.2f-%s.npy' % (freq, reim))
            if not os.path.exists(fpath):
                raise RuntimeError('missing file: %s' % fpath)
            log('Read %s' % fpath)
            if ksd_comm.rank == 0:
                rho_p = np.load(fpath)
            else:
                rho_p = None
            ksd.distribute_p(rho_p, rho_wrq[w, r])
            rho_wrq[w, r] /= kickstr

    omega_w = freq_w * eV_to_au
    for tag, mult_w in [('', np.ones_like(omega_w)),
                        ('-Iomega', -1j * omega_w),
                        ('-omega2', -1 * omega_w**2)]:
        # Do convolution time-by-time
        log('Start', tag)
        mult_wrq = mult_w[:, np.newaxis, np.newaxis]
        c = FourierConvolver(pulse, freq_w * eV_to_au, rho_wrq * mult_wrq)
        for t, time in enumerate(time_t):
            fpath = os.path.join(pulse_rho_dpath,
                                 't%09.1f%s.npy' % (time, tag))
            if os.path.exists(fpath):
                continue
            pulserho_rq = c.convolve(time, units='as')
            pulserho_q = pulserho_rq[0] + 1.0j * pulserho_rq[1]
            pulserho_p = ksd.collect_q(pulserho_q)
            log('Write %s' % fpath)
            if ksd_comm.rank == 0:
                np.save(fpath, pulserho_p)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, DirectoryPathType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=ExistingPathType)
    parser.add_argument('ksd_fpath', type=ExistingPathType)
    parser.add_argument('frho_dpath', type=ExistingPathType)
    parser.add_argument('pulse_fpath', type=ExistingPathType)
    parser.add_argument('pulse_rho_dpath', type=DirectoryPathType)
    parser.add_argument('--time', type=float, nargs='+')
    parser.add_argument('--maxtime', type=float, default=30e3)
    parser.add_argument('--freq_step', type=float, default=0.05)
    args = parser.parse_args()

    from pulse import read_pulse, convolution_freq_grid, pulse_time_grid
    pulse = read_pulse(args.pulse_fpath)
    freq_w = convolution_freq_grid(pulse, args.freq_step)
    time_t = pulse_time_grid(args.time, args.maxtime)

    pulse_rho(pulse, freq_w, time_t, args.gpw_fpath, args.ksd_fpath,
              args.frho_dpath, args.pulse_rho_dpath)
