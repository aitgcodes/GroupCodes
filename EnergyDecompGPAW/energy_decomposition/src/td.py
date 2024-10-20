import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))


def is_propagation_complete(dm_fpath, totaltime):
    from gpaw.tddft.units import au_to_as

    if os.path.exists(dm_fpath) and os.path.getsize(dm_fpath) > 0:
        data_ti = np.loadtxt(dm_fpath)
        finaltime = data_ti[-1, 0]
        if finaltime * au_to_as >= totaltime:
            return True
    return False


def td(gs_gpw_fpath,
       td_out_fpath,
       td_gpw_fpath,
       dm_fpath,
       fdm_fpath,
       kick,
       fxc,
       totaltime,
       fdm_freq_step):
    from ase.parallel import paropen
    from gpaw.mpi import world
    from gpaw.lcaotddft import LCAOTDDFT
    from gpaw.lcaotddft.dipolemomentwriter import DipoleMomentWriter
    from gpaw.lcaotddft.restartfilewriter import RestartFileWriter
    from gpaw.lcaotddft.densitymatrix import DensityMatrix
    from gpaw.lcaotddft.frequencydensitymatrix import FrequencyDensityMatrix
    from gpaw.tddft.units import au_to_as
    from gpaw.tddft.folding import frequencies
    from gpaw.utilities.timelimit import TimeLimiter, time_to_seconds

    from parallel_util import get_parallel

    # Basic settings
    kick_v = [0., 0., 0.]
    kick_v['xyz'.index(kick)] = 1e-5
    dt = 20.0
    totaltime -= 1e-6

    # Check if propagation is complete
    if is_propagation_complete(dm_fpath, totaltime):
        return

    # TimeLimiter settings from environment variables
    timelimit = time_to_seconds(os.environ['TIME_LIMIT'])
    timelimit -= time_to_seconds('40:0')
    timestart = float(os.environ['START_TIME'])

    # Calculate
    if not os.path.exists(td_gpw_fpath):
        # First time
        td_calc = LCAOTDDFT(gs_gpw_fpath, fxc=fxc,
                            parallel=get_parallel(world),
                            txt=paropen(td_out_fpath, 'w'))
        iterations = int(np.ceil(totaltime / dt))
        do_kick = True
        fdm_read_fpath = None
        fdm_freqs = frequencies(np.arange(0, 6 + 1e-3, fdm_freq_step),
                                None, None)
    else:
        # Restart
        td_calc = LCAOTDDFT(td_gpw_fpath,
                            parallel=get_parallel(world),
                            txt=paropen(td_out_fpath, 'a'))
        iterations = int(np.ceil((totaltime - td_calc.time * au_to_as) / dt))
        do_kick = False
        fdm_read_fpath = fdm_fpath
        fdm_freqs = None

    tl = TimeLimiter(td_calc, timelimit=timelimit, timestart=timestart)
    tl.reset('tddft')

    DipoleMomentWriter(td_calc, dm_fpath)
    RestartFileWriter(td_calc, td_gpw_fpath)
    dmat = DensityMatrix(td_calc)
    fdm = FrequencyDensityMatrix(td_calc, dmat, filename=fdm_read_fpath,
                                 frequencies=fdm_freqs)

    if do_kick:
        td_calc.absorption_kick(kick_v)

    td_calc.propagate(dt, iterations)
    td_calc.write(td_gpw_fpath, mode='all')
    fdm.write(fdm_fpath)


if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('gs_gpw_fpath', type=ExistingPathType)
    parser.add_argument('td_out_fpath', type=FilePathType)
    parser.add_argument('td_gpw_fpath', type=FilePathType)
    parser.add_argument('dm_fpath', type=FilePathType)
    parser.add_argument('fdm_fpath', type=FilePathType)
    parser.add_argument('--kick', default='x')
    parser.add_argument('--fxc', default=None)
    parser.add_argument('--totaltime', type=float, default=30e3)
    parser.add_argument('--fdm_freq_step', type=float, default=0.05)
    args = parser.parse_args()

    td(args.gs_gpw_fpath,
       args.td_out_fpath,
       args.td_gpw_fpath,
       args.dm_fpath,
       args.fdm_fpath,
       kick=args.kick,
       fxc=args.fxc,
       totaltime=args.totaltime,
       fdm_freq_step=args.fdm_freq_step)
