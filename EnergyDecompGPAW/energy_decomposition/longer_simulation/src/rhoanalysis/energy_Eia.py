import numpy as np
from rhoanalysis.base import BaseCalculator, build_filter


class EnergyCalculator(BaseCalculator):

    def run(self, flt=None):
        time_t = self.time_t

        flt_p = build_filter(self.ksd, flt)

        w_p = self.ksd.w_p[flt_p]

        ia_p = self.ksd.ia_p[flt_p]
        i_p = ia_p[:, 0]
        a_p = ia_p[:, 1]

        self.log('Calculate transition probability')
        imin = np.min(i_p)
        imax = np.max(i_p)
        amin = np.min(a_p)
        amax = np.max(a_p)

        Eia_t = np.zeros((len(time_t), imax - imin + 1, amax - amin + 1))

        weight_ia = np.zeros((imax - imin + 1, amax - amin + 1))

        read_keys = ['Q', 'P', 'dQ', 'dP', 'v']
        for t, Q_p, P_p, dQ_p, dP_p, v_p in self.read(read_keys, flt_p, v=0):
            PdQ = np.dot(P_p, dQ_p)
            QdP = np.dot(Q_p, dP_p)
            vQ = np.dot(v_p, Q_p)
            wQQ = np.dot(w_p * Q_p, Q_p)
            wPP = np.dot(w_p * P_p, P_p)

            weight_p = 0.5 * (PdQ - QdP - vQ)
            weight_ia[:,:] = 0.0
            np.add.at(weight_ia, (i_p - imin, a_p - amin), weight_p)
            Eia_t[t, :, :] = weight_ia

        return Eia_t

