import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def frho(gpw_fpath, ksd_fpath, fdm_fpath, frho_dpath):
    from gpaw.mpi import world
    from gpaw.mpi import SerialCommunicator
    from gpaw import GPAW
    from gpaw.lcaotddft.frequencydensitymatrix import \
        FrequencyDensityMatrixReader
    from gpaw.lcaotddft.ksdecomposition import KohnShamDecomposition
    from gpaw.tddft.units import au_to_eV
    from parallel_util import get_logger

    calc_comm = SerialCommunicator()
    loop_comm = world
    log = get_logger(calc_comm, loop_comm)

    calc = GPAW(gpw_fpath, txt=None, communicator=calc_comm)
    fdm = FrequencyDensityMatrixReader(fdm_fpath, calc.wfs.ksl, calc.wfs.kpt_u)
    ksd = KohnShamDecomposition(calc, ksd_fpath)

    log('Start')
    ffreq_w = fdm.freq_w
    reim_r = ['Re', 'Im']
    idx = -1
    for w, ffreq in enumerate(ffreq_w):
        for reim in reim_r:
            idx += 1
            if idx % loop_comm.size != loop_comm.rank:
                continue
            freq = ffreq.freq * au_to_eV
            folding = ffreq.folding.folding
            fname = 'w%05.2f-%s' % (freq, reim)
            if folding is not None:
                width = ffreq.folding.width * au_to_eV
                fname += '-%s-%.3f' % (folding, width)
            fname += '.npy'
            fpath = os.path.join(frho_dpath, fname)
            if os.path.exists(fpath):
                continue
            log('Calculate %s' % fname)
            rho_uMM = fdm.read_FDrho(reim, [w])[0]
            rho_p = ksd.transform(rho_uMM)[0]
            log('Write %s' % fpath)
            if calc_comm.rank == 0:
                np.save(fpath, rho_p)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, DirectoryPathType

    parser = argparse.ArgumentParser()
    parser.add_argument('gpw_fpath', type=ExistingPathType)
    parser.add_argument('ksd_fpath', type=ExistingPathType)
    parser.add_argument('fdm_fpath', type=ExistingPathType)
    parser.add_argument('frho_dpath', type=DirectoryPathType)
    args = parser.parse_args()

    frho(args.gpw_fpath, args.ksd_fpath, args.fdm_fpath, args.frho_dpath)
