import numpy as np
from gpaw.lcaotddft.ksdecomposition import gauss_ij
from rhoanalysis.base import BaseCalculator, build_filter


class WTECalculator(BaseCalculator):

    def run(self, weights_wfs, E_ia, outfpath, flt=None):
        time_t = self.time_t

        flt_p = build_filter(self.ksd, flt)
        eig_n, fermilevel = self.ksd.get_eig_n(zero_fermilevel=True)
        ia_p = self.ksd.ia_p[flt_p]
        i_p = ia_p[:, 0]
        a_p = ia_p[:, 1]

        imin = np.min(i_p)
        imax = np.max(i_p)
        amin = np.min(a_p)
        amax = np.max(a_p)

        # Weights of the occupied and unoccupied Kohn-Sham wave functions in the left 
        # and right regions of space
        wi_L = weights_wfs[0]
        wi_R = weights_wfs[1]
        wa_L = weights_wfs[2]
        wa_R = weights_wfs[3]

        # Weighted transition energy
        E_LL_t = np.zeros(len(time_t))
        E_RR_t = np.zeros(len(time_t))
        E_LR_t = np.zeros(len(time_t))
        E_RL_t = np.zeros(len(time_t))

        
        for t in range(len(time_t)):
            for i in range(imin, imax+1):
                for a in range(amin, amax+1):
                    E_LL_t[t] += wi_L[i-imin]*wa_L[a-amin]*E_ia[t,i-imin,a-amin]
                    E_RR_t[t] += wi_R[i-imin]*wa_R[a-amin]*E_ia[t,i-imin,a-amin]
                    E_LR_t[t] += wi_L[i-imin]*wa_R[a-amin]*E_ia[t,i-imin,a-amin]
                    E_RL_t[t] += wi_R[i-imin]*wa_L[a-amin]*E_ia[t,i-imin,a-amin]

        if self.world.rank == 0:
            np.savez_compressed(outfpath,time_t=time_t,E_LL_t=E_LL_t,E_RR_t=E_RR_t,E_LR_t=E_LR_t,E_RL_t=E_RL_t )
