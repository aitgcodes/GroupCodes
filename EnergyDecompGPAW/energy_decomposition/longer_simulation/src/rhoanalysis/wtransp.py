import numpy as np
from gpaw.lcaotddft.ksdecomposition import gauss_ij
from rhoanalysis.base import BaseCalculator, build_filter


class WTPCalculator(BaseCalculator):

    def run(self, energy_o, energy_u, sigma, weights_wfs, P_ia, outfpath, flt=None):
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

        # Total probability of a partial process (in %)
        N_LR_t = np.zeros(len(time_t))
        N_LL_t = np.zeros(len(time_t))
        N_RR_t = np.zeros(len(time_t))
        N_RL_t = np.zeros(len(time_t))

        P_LR_t = np.zeros(len(time_t))
        P_LL_t = np.zeros(len(time_t))
        P_RR_t = np.zeros(len(time_t))
        P_RL_t = np.zeros(len(time_t))


        N_t = np.zeros(len(time_t))

        
        for t in range(len(time_t)):
            for i in range(imin, imax+1):
                for a in range(amin, amax+1):
                    N_t[t] += P_ia[t,i-imin,a-amin]
                    N_LR_t[t] += wi_L[i-imin]*wa_R[a-amin]*P_ia[t,i-imin,a-amin]
                    N_LL_t[t] += wi_L[i-imin]*wa_L[a-amin]*P_ia[t,i-imin,a-amin]
                    N_RR_t[t] += wi_R[i-imin]*wa_R[a-amin]*P_ia[t,i-imin,a-amin]
                    N_RL_t[t] += wi_R[i-imin]*wa_L[a-amin]*P_ia[t,i-imin,a-amin]

            P_LR_t[t] = (N_LR_t[t]/N_t[t])*100
            P_LL_t[t] = (N_LL_t[t]/N_t[t])*100
            P_RR_t[t] = (N_RR_t[t]/N_t[t])*100
            P_RL_t[t] = (N_RL_t[t]/N_t[t])*100



        if self.world.rank == 0:
            np.savez_compressed(outfpath, time_t=time_t,N_t=N_t,N_LR_t=N_LR_t,N_LL_t=N_LL_t,
                    N_RR_t=N_RR_t,N_RL_t=N_RL_t,
                    P_LR_t=P_LR_t,P_LL_t=P_LL_t,P_RR_t=P_RR_t,P_RL_t=P_RL_t)
