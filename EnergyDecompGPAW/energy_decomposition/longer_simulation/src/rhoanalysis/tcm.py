import os
import numpy as np
from gpaw.lcaotddft.ksdecomposition import gauss_ij
from rhoanalysis.base import BaseCalculator, get_weight_p, get_keys


class TCMCalculator(BaseCalculator):

    def run(self, energy_o, energy_u, sigma, outdpath, weight):
        eig_n, fermilevel = self.ksd.get_eig_n(zero_fermilevel=True)
        flt_p = self.ksd.filter_by_x_ia(eig_n, energy_o, energy_u, 8 * sigma)

        self.log('Calculate gauss_ij')
        G_po = gauss_ij(eig_n[self.ksd.ia_p[flt_p, 0]], energy_o, sigma)
        G_pu = gauss_ij(eig_n[self.ksd.ia_p[flt_p, 1]], energy_u, sigma)

        w_p = self.ksd.w_p[flt_p]

        read_keys = get_keys(weight)
        for t, data_k in self.read(read_keys, flt_p, v=0, yield_dict=True):
            outfpath = os.path.join(outdpath, 't%09.1f.npz' % self.time_t[t])
            weight_p = get_weight_p(weight, w_p=w_p, **data_k)
            tcm_ou = np.dot(G_po.T * weight_p, G_pu)
            if self.calc_comm.rank == 0:
                np.savez_compressed(outfpath, time=self.time_t[t],
                                    energy_o=energy_o, energy_u=energy_u,
                                    sigma=sigma, tcm_ou=tcm_ou)
