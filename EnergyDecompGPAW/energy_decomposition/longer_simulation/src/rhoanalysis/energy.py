import numpy as np
from rhoanalysis.base import BaseCalculator, build_filter


class EnergyCalculator(BaseCalculator):

    def run(self, outfpath, flt=None):
        time_t = self.time_t
        flt_p = build_filter(self.ksd, flt)

        w_p = self.ksd.w_p[flt_p]

        E_t = np.zeros(len(time_t))
        Ec_t = np.zeros(len(time_t))
        Eq_t = np.zeros(len(time_t))
        Ep_t = np.zeros(len(time_t))

        read_keys = ['Q', 'P', 'dQ', 'dP', 'v']
        for t, Q_p, P_p, dQ_p, dP_p, v_p in self.read(read_keys, flt_p, v=0):
            PdQ = np.dot(P_p, dQ_p)
            QdP = np.dot(Q_p, dP_p)
            vQ = np.dot(v_p, Q_p)
            wQQ = np.dot(w_p * Q_p, Q_p)
            wPP = np.dot(w_p * P_p, P_p)
            E_t[t] = 0.5 * (PdQ - QdP - vQ)
            Ec_t[t] = -0.5 * (wQQ + QdP + vQ)
            Ep_t[t] = 0.5 * wPP
            Eq_t[t] = E_t[t] - Ep_t[t]

        self.loop_comm.sum(E_t, 0)
        self.loop_comm.sum(Ec_t, 0)
        self.loop_comm.sum(Eq_t, 0)
        self.loop_comm.sum(Ep_t, 0)

        if self.world.rank == 0:
            np.savez_compressed(outfpath, time_t=time_t,
                                E_t=E_t, Ec_t=Ec_t, Eq_t=Eq_t, Ep_t=Ep_t)
