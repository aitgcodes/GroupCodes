import numpy as np
from rhoanalysis.base import BaseCalculator, build_filter


class HCAtomAmountCalculator(BaseCalculator):

    def __init__(self, *args, **kwargs):
        self.weight_am = kwargs.pop('weight_am')
        self.n_m = kwargs.pop('n_m')
        self.d_m = kwargs.pop('d_m')
        self.deg_n = kwargs.pop('deg_n')
        BaseCalculator.__init__(self, *args, **kwargs)

    def run(self, outfpath, flt=None):
        time_t = self.time_t
        Na = self.weight_am.shape[0]
        Nn = len(self.deg_n)
        fi_ta = np.zeros((len(time_t), Na))
        fa_ta = np.zeros((len(time_t), Na))

        flt_p = build_filter(self.ksd, flt)
        ia_p = self.ksd.ia_p[flt_p]
        i_p = ia_p[:, 0]
        a_p = ia_p[:, 1]
        imin = np.min(i_p)
        imax = np.max(i_p)
        amin = np.min(a_p)
        amax = np.max(a_p)

        def M_ia_from_M_p(M_p):
            M_ia = np.zeros((imax - imin + 1, amax - amin + 1),
                            dtype=M_p.dtype)
            for M, (i, a) in zip(M_p, ia_p):
                M_ia[i - imin, a - amin] = M
            return M_ia

        # Build weight matrix
        Nd = np.max(self.deg_n)
        weight_dan = []
        for d in range(Nd):
            weight_an = np.zeros((Na, Nn))
            flt_m = self.d_m == d
            weight_an[:, self.n_m[flt_m]] = self.weight_am[:, flt_m]
            weight_dan.append(weight_an)

        read_keys = ['Q', 'P']
        for t, Q_p, P_p in self.read(read_keys, flt_p, v=0):
            Q_ia = M_ia_from_M_p(Q_p)
            P_ia = M_ia_from_M_p(P_p)

            # Holes
            M_ii = 0.5 * (np.dot(Q_ia, Q_ia.T) + np.dot(P_ia, P_ia.T))
            # Loop over diagonals
            for d, weight_an in enumerate(weight_dan):
                # Multiplier is 1.0 for the main diagonal and
                # 2.0 for diagonals above the main diagonal
                mult = 1.0 + float(d != 0)
                weight_ai = weight_an[:, (imin + d):(imax + 1)]
                if np.any(np.isnan(weight_an[:, (imin + d):(imax + 1)])):
                    self.log('not enough atomweights calculated for holes')
                fi_ta[t, :] += mult * np.dot(weight_ai, np.diag(M_ii, d))

            # Electrons
            M_aa = 0.5 * (np.dot(Q_ia.T, Q_ia) + np.dot(P_ia.T, P_ia))
            # Loop over diagonals
            for d, weight_an in enumerate(weight_dan):
                # Multiplier is 1.0 for the main diagonal and
                # 2.0 for diagonals above the main diagonal
                mult = 1.0 + float(d != 0)
                weight_aa = weight_an[:, (amin + d):(amax + 1)]
                if np.any(np.isnan(weight_an[:, (amin + d):(amax + 1)])):
                    self.log('not enough atomweights calculated for electrons')
                fa_ta[t, :] += mult * np.dot(weight_aa, np.diag(M_aa, d))

        self.loop_comm.sum(fi_ta, 0)
        self.loop_comm.sum(fa_ta, 0)

        if self.world.rank == 0:
            np.savez_compressed(outfpath, time_t=time_t,
                                fi_ta=fi_ta, fa_ta=fa_ta)
