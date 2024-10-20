import numpy as np
from gpaw.lcaotddft.ksdecomposition import gauss_ij
from rhoanalysis.base import BaseCalculator, build_filter


class WTPDistCalculator(BaseCalculator):

    def run(self, energy_o, energy_u, sigma, weights_wfs, P_ia, outfpath, flt=None):
        time_t = self.time_t

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

        #print("shape of G_io =", G_io.shape)
        #print("shape of G_au =", G_au.shape)


        # Weights of the occupied and unoccupied Kohn-Sham wave functions in the left 
        # and right regions of space
        wi_L = weights_wfs[0]
        wi_R = weights_wfs[1]
        wa_L = weights_wfs[2]
        wa_R = weights_wfs[3]

        Pe_LR_tu = np.zeros((len(time_t), len(energy_u)))
        Pe_LL_tu = np.zeros((len(time_t), len(energy_u)))
        Pe_RR_tu = np.zeros((len(time_t), len(energy_u)))
        Pe_RL_tu = np.zeros((len(time_t), len(energy_u)))
        
        Ph_LR_to = np.zeros((len(time_t), len(energy_o)))
        Ph_LL_to = np.zeros((len(time_t), len(energy_o)))
        Ph_RR_to = np.zeros((len(time_t), len(energy_o)))
        Ph_RL_to = np.zeros((len(time_t), len(energy_o)))

        # Weighted transition probability of electrons
        for t in range(len(time_t)):
            for eu in range(len(energy_u)):
                for i in range(imin, imax+1):
                    for a in range(amin, amax+1):
                        Pe_LR_tu[t, eu] += wi_L[i-imin]*wa_R[a-amin]*P_ia[t,i-imin,a-amin]*G_au[a-amin,eu]
                        Pe_LL_tu[t, eu] += wi_L[i-imin]*wa_L[a-amin]*P_ia[t,i-imin,a-amin]*G_au[a-amin,eu]
                        Pe_RR_tu[t, eu] += wi_R[i-imin]*wa_R[a-amin]*P_ia[t,i-imin,a-amin]*G_au[a-amin,eu]
                        Pe_RL_tu[t, eu] += wi_R[i-imin]*wa_L[a-amin]*P_ia[t,i-imin,a-amin]*G_au[a-amin,eu]
        # Weighted transition probability of holes
        for t in range(len(time_t)):
            for eo in range(len(energy_o)):
                for i in range(imin, imax+1):
                    for a in range(amin, amax+1):
                        Ph_LR_to[t, eo] += wi_L[i-imin]*wa_R[a-amin]*P_ia[t,i-imin,a-amin]*G_io[i-imin,eo]
                        Ph_LL_to[t, eo] += wi_L[i-imin]*wa_L[a-amin]*P_ia[t,i-imin,a-amin]*G_io[i-imin,eo]
                        Ph_RR_to[t, eo] += wi_R[i-imin]*wa_R[a-amin]*P_ia[t,i-imin,a-amin]*G_io[i-imin,eo]
                        Ph_RL_to[t, eo] += wi_R[i-imin]*wa_L[a-amin]*P_ia[t,i-imin,a-amin]*G_io[i-imin,eo]

        if self.world.rank == 0:
            np.savez_compressed(outfpath, time_t=time_t,
                    energy_o=energy_o,Ph_LR_to=Ph_LR_to,Ph_LL_to=Ph_LL_to,
                    Ph_RR_to=Ph_RR_to,Ph_RL_to=Ph_RL_to,
                    energy_u=energy_u, Pe_LR_tu=Pe_LR_tu,Pe_LL_tu=Pe_LL_tu,
                    Pe_RR_tu=Pe_RR_tu,Pe_RL_tu=Pe_RL_tu)
