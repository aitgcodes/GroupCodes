import numpy as np
from gpaw.lcaotddft.ksdecomposition import gauss_ij
from gpaw.tddft.units import au_to_eV
from rhoanalysis.base import BaseCalculator, get_weight_p, get_keys


class OmegaDistCalculator(BaseCalculator):

    def run(self, energy_e, sigma, outfpath, weight):
        time_t = self.time_t
        dist_te = np.zeros((len(time_t), len(energy_e)))

        flt_p = self.ksd.filter_by_x_p(self.ksd.w_p * au_to_eV,
                                       energy_e, 8 * sigma)
        w_p = self.ksd.w_p[flt_p]

        self.log('Calculate gauss_ij')
        G_pe = gauss_ij(w_p * au_to_eV, energy_e, sigma)

        read_keys = get_keys(weight)
        for t, data_k in self.read(read_keys, flt_p, v=0, yield_dict=True):
            weight_p = get_weight_p(weight, w_p=w_p, **data_k)
            dist_e = np.dot(G_pe.T, weight_p)
            dist_te[t, :] = dist_e

        self.loop_comm.sum(dist_te, 0)

        if self.world.rank == 0:
            np.savez_compressed(outfpath, time_t=time_t,
                                energy_e=energy_e, dist_te=dist_te)
