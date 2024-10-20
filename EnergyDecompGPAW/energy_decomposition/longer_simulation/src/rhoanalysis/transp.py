import numpy as np
from gpaw.lcaotddft.ksdecomposition import gauss_ij
from rhoanalysis.base import BaseCalculator, build_filter


class TransPCalculator(BaseCalculator):

    #def run(self, outfpath, flt=None):
    def run(self, flt=None):
        time_t = self.time_t

        flt_p = build_filter(self.ksd, flt)
        eig_n, fermilevel = self.ksd.get_eig_n(zero_fermilevel=True)
        ia_p = self.ksd.ia_p[flt_p]
        i_p = ia_p[:, 0]
        a_p = ia_p[:, 1]

        self.log('Calculate transition probability')
        imin = np.min(i_p)
        imax = np.max(i_p)
        amin = np.min(a_p)
        amax = np.max(a_p)

        print("imin = ", imin)
        print("imax = ", imax)
        print("amin = ", amin)
        print("amax = ", amax)

        transp = np.zeros((len(time_t), imax - imin + 1, amax - amin + 1))
        
        weight_ia = np.zeros((imax - imin + 1, amax - amin + 1))


        read_keys = ['Q', 'P']
        for t, Q_p, P_p in self.read(read_keys, flt_p, v=0):
            weight_p = 0.5 * (Q_p**2 + P_p**2)
            weight_ia[:,:] = 0.0
            np.add.at(weight_ia, (i_p - imin, a_p - amin), weight_p)
            transp[t, :, :] = weight_ia


        #if self.world.rank == 0:
        #    np.savez_compressed(outfpath, time_t=time_t,
        #                        transp=transp)

        return transp


