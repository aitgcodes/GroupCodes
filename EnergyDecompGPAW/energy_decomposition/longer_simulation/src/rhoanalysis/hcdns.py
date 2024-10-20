import numpy as np
import ase.io
from gpaw.lcaotddft.densitymatrix import get_density
from rhoanalysis.base import BaseCalculator, build_filter


class HCDensityCalculator(BaseCalculator):

    def __init__(self, *args, **kwargs):
        self.deg_n = kwargs.pop('deg_n')
        BaseCalculator.__init__(self, *args, **kwargs)
        # Additional initialization in order to calculate density
        self.calc.initialize_positions()

    def run(self, elho, outfpath, flt=None):
        time_t = self.time_t
        if len(time_t) != 1:
            raise RuntimeError('output file is only for one time instance')

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

        # Choose electron or holes
        if elho == 'holes':
            jmin = imin
            jmax = imax
        elif elho == 'electrons':
            jmin = amin
            jmax = amax

        # Build filter for density matrix
        flt_jj = np.zeros((jmax - jmin + 1, jmax - jmin + 1), dtype=int)
        deg_j = self.deg_n[jmin:(jmax + 1)]
        if deg_j[0] != 1:
            raise RuntimeError('adjust filter limit to start from '
                               'a non-degenerate state')
        for j, deg in enumerate(deg_j):
            flt_jj[j, j] = 1
            for k in range(1, deg):
                flt_jj[j - k, j] = 1
                flt_jj[j, j - k] = 1

        # Wave functions
        C0_jM = self.ksd.C0_unM[0][jmin:(jmax + 1)]

        read_keys = ['Q', 'P']
        for t, Q_p, P_p in self.read(read_keys, flt_p, v=0):
            Q_ia = M_ia_from_M_p(Q_p)
            P_ia = M_ia_from_M_p(P_p)

            if elho == 'holes':
                M_jj = 0.5 * (np.dot(Q_ia, Q_ia.T) + np.dot(P_ia, P_ia.T))
            elif elho == 'electrons':
                M_jj = 0.5 * (np.dot(Q_ia.T, Q_ia) + np.dot(P_ia.T, P_ia))
            M_jj = M_jj * flt_jj

            # Construct density
            rho_MM = np.dot(C0_jM.T, np.dot(M_jj, C0_jM.conj()))
            rho_MM = 0.5 * (rho_MM + rho_MM.T)
            rho_g = get_density(rho_MM, self.calc.wfs, self.ksd.density,
                                'comp', 0)
            gd = self.ksd.density.finegd
            total = np.trace(M_jj)
            rerr = np.absolute(gd.integrate(rho_g) - total) / total
            self.log('rerr: %e' % rerr)

            big_g = gd.collect(rho_g)
            if self.calc_comm.rank == 0:
                pad_g = np.zeros(np.array(big_g.shape) + 1)
                pad_g[1:, 1:, 1:] = big_g
                ase.io.write(outfpath, self.calc.atoms, data=pad_g)
