import numpy as np
from gpaw.lcaotddft.ksdecomposition import gauss_ij
from rhoanalysis.base import BaseCalculator


class HCAtomDistCalculator(BaseCalculator):

    def __init__(self, *args, **kwargs):
        self.weight_am = kwargs.pop('weight_am')
        self.n_m = kwargs.pop('n_m')
        self.d_m = kwargs.pop('d_m')
        self.deg_n = kwargs.pop('deg_n')
        BaseCalculator.__init__(self, *args, **kwargs)

    def run(self, energy_o, energy_u, sigma, outfpath, aproj_i):
        time_t = self.time_t
        dist_to = np.zeros((len(time_t), len(energy_o)))
        dist_tu = np.zeros((len(time_t), len(energy_u)))

        eig_n, fermilevel = self.ksd.get_eig_n(zero_fermilevel=True)

        # Build weight matrix
        Nd = np.max(self.deg_n)
        weight_m = np.sum(self.weight_am[aproj_i], axis=0)
        flt_n = np.logical_and(eig_n >= energy_o.min() - 8 * sigma,
                               eig_n <= energy_u.max() + 8 * sigma)
        weight_dn = []
        for d in range(Nd):
            weight_n = np.zeros_like(eig_n)
            flt_m = self.d_m == d
            weight_n[self.n_m[flt_m]] = weight_m[flt_m]
            if np.any(np.isnan(weight_n[flt_n])):
                raise RuntimeError('not enough atom weights calculated')
            weight_n = np.nan_to_num(weight_n)
            weight_dn.append(weight_n)

        self.log('Calculate gauss_ij')
        imin, imax, amin, amax = self.ksd.ialims()
        G_io = gauss_ij(eig_n[imin:(imax + 1)], energy_o, sigma)
        G_au = gauss_ij(eig_n[amin:(amax + 1)], energy_u, sigma)

        read_keys = ['Q', 'P']
        for t, Q_p, P_p in self.read(read_keys, v=0):
            Q_ia = self.ksd.M_ia_from_M_p(Q_p)
            P_ia = self.ksd.M_ia_from_M_p(P_p)

            # Holes
            M_ii = 0.5 * (np.dot(Q_ia, Q_ia.T) + np.dot(P_ia, P_ia.T))
            # Loop over diagonals
            for d, weight_n in enumerate(weight_dn):
                # Multiplier is 1.0 for the main diagonal and
                # 2.0 for diagonals above the main diagonal
                mult = 1.0 + float(d != 0)
                weight_i = weight_n[(imin + d):(imax + 1)] * np.diag(M_ii, d)
                dist_to[t, :] += mult * np.dot(weight_i, G_io[d:])

            # Electrons
            M_aa = 0.5 * (np.dot(Q_ia.T, Q_ia) + np.dot(P_ia.T, P_ia))
            # Loop over diagonals
            for d, weight_n in enumerate(weight_dn):
                # Multiplier is 1.0 for the main diagonal and
                # 2.0 for diagonals above the main diagonal
                mult = 1.0 + float(d != 0)
                weight_a = weight_n[(amin + d):(amax + 1)] * np.diag(M_aa, d)
                dist_tu[t, :] += mult * np.dot(weight_a, G_au[d:])

        self.loop_comm.sum(dist_to, 0)
        self.loop_comm.sum(dist_tu, 0)

        if self.world.rank == 0:
            np.savez_compressed(outfpath, time_t=time_t,
                                energy_o=energy_o, dist_to=dist_to,
                                energy_u=energy_u, dist_tu=dist_tu)
