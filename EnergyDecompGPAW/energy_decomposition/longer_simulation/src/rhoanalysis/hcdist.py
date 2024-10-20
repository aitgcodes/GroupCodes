import numpy as np
from gpaw.lcaotddft.ksdecomposition import gauss_ij
from rhoanalysis.base import BaseCalculator, build_filter


class HCDistCalculator(BaseCalculator):

    def run(self, energy_o, energy_u, sigma, outfpath, flt=None):
        time_t = self.time_t
        dist_to = np.zeros((len(time_t), len(energy_o)))
        dist_tu = np.zeros((len(time_t), len(energy_u)))

        flt_p = build_filter(self.ksd, flt)
        eig_n, fermilevel = self.ksd.get_eig_n(zero_fermilevel=True)
        ia_p = self.ksd.ia_p[flt_p]
        i_p = ia_p[:, 0]
        a_p = ia_p[:, 1]

        self.log('Calculate gauss_ij')
        imin = np.min(i_p)
        imax = np.max(i_p)
        amin = np.min(a_p)
        amax = np.max(a_p)
        G_io = gauss_ij(eig_n[imin:(imax + 1)], energy_o, sigma)
        G_au = gauss_ij(eig_n[amin:(amax + 1)], energy_u, sigma)

        weight_i = np.zeros(imax - imin + 1)
        weight_a = np.zeros(amax - amin + 1)

        read_keys = ['Q', 'P']
        for t, Q_p, P_p in self.read(read_keys, flt_p, v=0):
            weight_p = 0.5 * (Q_p**2 + P_p**2)
            weight_i[:] = 0.0
            weight_a[:] = 0.0
            np.add.at(weight_i, i_p - imin, weight_p)
            np.add.at(weight_a, a_p - amin, weight_p)
            dist_o = np.dot(weight_i, G_io)
            dist_u = np.dot(weight_a, G_au)
            dist_to[t, :] = dist_o
            dist_tu[t, :] = dist_u

        self.loop_comm.sum(dist_to, 0)
        self.loop_comm.sum(dist_tu, 0)

        if self.world.rank == 0:
            np.savez_compressed(outfpath, time_t=time_t,
                                energy_o=energy_o, dist_to=dist_to,
                                energy_u=energy_u, dist_tu=dist_tu)
