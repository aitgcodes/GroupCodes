import os
import sys
import numpy as np

from gpaw.mpi import world
from gpaw.mpi import SerialCommunicator
from gpaw import GPAW
from gpaw.lcaotddft.ksdecomposition import KohnShamDecomposition
from gpaw.tddft.units import as_to_au, au_to_eV


def build_single_filter(ksd, key, lims):
    eig_n, fermilevel = ksd.get_eig_n(zero_fermilevel=True)

    if key in ['e', 'u', 'a']:
        x_p = eig_n[ksd.ia_p[:, 1]]
    elif key in ['h', 'o', 'i']:
        x_p = eig_n[ksd.ia_p[:, 0]]
    elif key in ['w']:
        x_p = ksd.w_p * au_to_eV
    else:
        raise RuntimeError('unknown filter: {}'.format(key))
    return np.logical_and(x_p >= lims[0], x_p < lims[1])


def build_filter(ksd, flt_list):
    if flt_list is None:
        flt_p = slice(None)
    else:
        key, lims = flt_list.pop()
        flt_p = build_single_filter(ksd, key, lims)
        while len(flt_list) > 0:
            op = flt_list.pop()
            key, lims = flt_list.pop()
            flt1_p = build_single_filter(ksd, key, lims)
            if op == 'or':
                flt_p = np.logical_or(flt_p, flt1_p)
            elif op == 'and':
                flt_p = np.logical_and(flt_p, flt1_p)
            else:
                raise RuntimeError('unknown operator: {}'.format(op))
    return flt_p


def get_keys(name):
    if name in ['occupation']:
        keys = ['Q', 'P']
    elif name in ['energy']:
        keys = ['Q', 'P', 'dQ', 'dP', 'v']
    elif name in ['coulombenergy']:
        keys = ['Q', 'dP', 'v']
    elif name in ['rate_absorption']:
        keys = ['dQ', 'v']
    elif name in ['rate_interaction']:
        keys = ['Q', 'P', 'dQ', 'ddQ', 'ddP', 'v', 'dv']
    else:
        raise RuntimeError('Unknown weight: {}'.format(name))
    return keys


def get_weight_p(name, w_p=None,
                 Q_p=None, P_p=None, v_p=None,
                 dQ_p=None, dP_p=None, dv_p=None,
                 ddQ_p=None, ddP_p=None):
    if name in ['occupation']:
        weight_p = P_p**2 + Q_p**2
    elif name in ['energy']:
        weight_p = 0.5 * (P_p * dQ_p - Q_p * dP_p - v_p * Q_p)
    elif name in ['coulombenergy']:
        weight_p = -0.5 * (w_p * Q_p**2 + Q_p * dP_p + v_p * Q_p)
    elif name in ['rate_absorption']:
        weight_p = -v_p * dQ_p
    elif name in ['rate_interaction']:
        weight_p = 0.5 * (P_p * ddQ_p - Q_p * ddP_p + v_p * dQ_p - dv_p * Q_p)
    else:
        raise RuntimeError('Unknown weight: {}'.format(name))
    return weight_p


class BaseCalculator(object):

    def __init__(self, time_t, pulse, gpw_fpath, ksd_fpath, pulse_rho_dpath,
                 continue_on_fail=False):
        self.calc_comm = SerialCommunicator()
        self.loop_comm = world
        self.world = world
        self.rho_dpath = pulse_rho_dpath
        self.time_t = time_t
        self.autime_t = time_t * as_to_au
        self.pulse = pulse
        self.continue_on_fail = continue_on_fail

        # Load ksd
        self.calc = GPAW(gpw_fpath, txt=None, communicator=self.calc_comm)
        self.ksd = KohnShamDecomposition(self.calc, ksd_fpath)

    def read(self, keys=['Q', 'P', 'dQ', 'dP', 'ddQ', 'ddP', 'v', 'dv'],
             flt_p=slice(None), v=0, yield_dict=False):
        f_p = self.ksd.f_p[flt_p]
        if 'v' in keys or 'dv' in keys:
            v0_p = self.ksd.dm_vp[v, flt_p] * np.sqrt(2 * f_p)

        tag_s = ['', '-Iomega', '-omega2']
        t_i = range(self.loop_comm.rank, len(self.time_t), self.loop_comm.size)
        for t in t_i:
            time = self.time_t[t]
            self.log('%10.1f as' % time)
            data_k = {}

            if any([k in keys for k in ['Q', 'P']]):
                rho_p = self.read_rho(time, tag_s[0], flt_p)
                if rho_p is None:
                    continue
                if 'Q' in keys:
                    Q_p = rho_p.real / np.sqrt(2 * f_p)
                    data_k['Q_p'] = Q_p
                if 'P' in keys:
                    P_p = rho_p.imag / np.sqrt(2 * f_p)
                    data_k['P_p'] = P_p
                rho_p = None

            if any([k in keys for k in ['dQ', 'dP']]):
                drho_p = self.read_rho(time, tag_s[1], flt_p)
                if drho_p is None:
                    continue
                if 'dQ' in keys:
                    dQ_p = drho_p.real / np.sqrt(2 * f_p)
                    data_k['dQ_p'] = dQ_p
                if 'dP' in keys:
                    dP_p = drho_p.imag / np.sqrt(2 * f_p)
                    data_k['dP_p'] = dP_p
                drho_p = None

            if any([k in keys for k in ['ddQ', 'ddP']]):
                ddrho_p = self.read_rho(time, tag_s[2], flt_p)
                if ddrho_p is None:
                    continue
                if 'ddQ' in keys:
                    ddQ_p = ddrho_p.real / np.sqrt(2 * f_p)
                    data_k['ddQ_p'] = ddQ_p
                if 'ddP' in keys:
                    ddP_p = ddrho_p.imag / np.sqrt(2 * f_p)
                    data_k['ddP_p'] = ddP_p
                ddrho_p = None

            if 'v' in keys:
                v_p = v0_p * self.pulse.strength(self.autime_t[t])
                data_k['v_p'] = v_p

            if 'dv' in keys:
                dv_p = v0_p * self.pulse.derivative(self.autime_t[t])
                data_k['dv_p'] = dv_p

            if yield_dict:
                yield t, data_k
            else:
                yield [t] + [data_k[k + '_p'] for k in keys]

    def read_rho(self, time, tag, flt_p):
        fpath = os.path.join(self.rho_dpath, 't%09.1f%s.npy' % (time, tag))
        if not os.path.exists(fpath):
            if self.continue_on_fail:
                return None
            else:
                raise RuntimeError('File missing: %s' % fpath)
        rho_p = np.load(fpath)[flt_p]
        return rho_p

    def run(self):
        raise NotImplementedError()

    def log(self, msg):
        r = self.loop_comm.rank
        s = self.loop_comm.size
        print('[%04d/%04d] %s' % (r, s, msg))
        sys.stdout.flush()
