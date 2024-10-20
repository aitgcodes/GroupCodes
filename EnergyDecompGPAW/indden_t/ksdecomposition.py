import numpy as np

from ase.units import Hartree, Bohr

from ase.io.ulm import Reader
from gpaw.io import Writer
from gpaw.external import ConstantElectricField
from gpaw.lcaotddft.hamiltonian import KickHamiltonian
from gpaw.lcaotddft.utilities import collect_MM
from gpaw.lcaotddft.utilities import distribute_nM
from gpaw.lcaotddft.utilities import read_uMM
from gpaw.lcaotddft.utilities import write_uMM
from gpaw.lcaotddft.utilities import read_uX, write_uX
from gpaw.utilities.scalapack import \
    pblas_simple_gemm, pblas_simple_hemm, scalapack_tri2full
from gpaw.utilities.tools import tri2full


def gauss_ij(energy_i, energy_j, sigma):
    denergy_ij = energy_i[:, np.newaxis] - energy_j[np.newaxis, :]
    norm = 1.0 / (sigma * np.sqrt(2 * np.pi))
    return norm * np.exp(-0.5 * denergy_ij**2 / sigma**2)


def get_bfs_maps(calc):
    # Construct maps
    # a_M: M -> atom index a
    # l_M: M -> angular momentum l
    a_M = []
    l_M = []
    M = 0
    for a, sphere in enumerate(calc.wfs.basis_functions.sphere_a):
        for j, spline in enumerate(sphere.spline_j):
            l = spline.get_angular_momentum_number()
            for _ in range(2 * l + 1):
                a_M.append(a)
                l_M.append(l)
                M += 1
    a_M = np.array(a_M)
    l_M = np.array(l_M)
    return a_M, l_M


class KohnShamDecomposition:
    version = 1
    ulmtag = 'KSD'
    readwrite_attrs = ['fermilevel', 'only_ia', 'w_p', 'f_p', 'ia_p',
                       'P_p', 'dm_vp', 'a_M', 'l_M']

    def __init__(self, paw=None, filename=None):
        self.filename = filename
        self.has_initialized = False
        self.reader = None
        if paw is not None:
            self.world = paw.world
            self.log = paw.log
            self.ksl = paw.wfs.ksl
            self.kd = paw.wfs.kd
            self.bd = paw.wfs.bd
            self.kpt_u = paw.wfs.kpt_u
            self.density = paw.density
            self.comm = paw.comms['K']

            if len(paw.wfs.kpt_u) > 1:
                raise RuntimeError('K-points are not fully supported')

        if filename is not None:
            self.read(filename)
            return

    def initialize(self, paw, min_occdiff=1e-3, only_ia=True):
        if self.has_initialized:
            return
        paw.initialize_positions()
        # paw.set_positions()

        assert self.bd.nbands == self.ksl.nao
        self.only_ia = only_ia

        if not self.ksl.using_blacs and self.bd.comm.size > 1:
            raise RuntimeError('Band parallelization without scalapack '
                               'is not supported')

        if self.kd.gamma:
            self.C0_dtype = float
        else:
            self.C0_dtype = complex

        # Take quantities
        self.fermilevel = paw.wfs.fermi_level
        self.S_uMM = []
        self.C0_unM = []
        self.eig_un = []
        self.occ_un = []
        for kpt in paw.wfs.kpt_u:
            S_MM = kpt.S_MM
            assert np.max(np.absolute(S_MM.imag)) == 0.0
            S_MM = np.ascontiguousarray(S_MM.real)
            if self.ksl.using_blacs:
                scalapack_tri2full(self.ksl.mmdescriptor, S_MM)
            self.S_uMM.append(S_MM)

            C_nM = kpt.C_nM
            if self.C0_dtype == float:
                assert np.max(np.absolute(C_nM.imag)) == 0.0
                C_nM = np.ascontiguousarray(C_nM.real)
            C_nM = distribute_nM(self.ksl, C_nM)
            self.C0_unM.append(C_nM)

            eig_n = paw.wfs.collect_eigenvalues(kpt.k, kpt.s)
            occ_n = paw.wfs.collect_occupations(kpt.k, kpt.s)
            self.eig_un.append(eig_n)
            self.occ_un.append(occ_n)

        self.a_M, self.l_M = get_bfs_maps(paw)
        self.atoms = paw.atoms

        # TODO: do the rest of the function with K-points

        # Construct p = (i, a) pairs
        u = 0
        eig_n = self.eig_un[u]
        occ_n = self.occ_un[u]
        C0_nM = self.C0_unM[u]

        if self.comm.rank == 0:
            Nn = self.bd.nbands

            f_p = []
            w_p = []
            i_p = []
            a_p = []
            ia_p = []
            i0 = 0
            for i in range(i0, Nn):
                if only_ia:
                    a0 = i + 1
                else:
                    a0 = 0
                for a in range(a0, Nn):
                    f = occ_n[i] - occ_n[a]
                    if only_ia and f < min_occdiff:
                        continue
                    w = eig_n[a] - eig_n[i]
                    f_p.append(f)
                    w_p.append(w)
                    i_p.append(i)
                    a_p.append(a)
                    ia_p.append((i, a))
            f_p = np.array(f_p)
            w_p = np.array(w_p)
            i_p = np.array(i_p, dtype=int)
            a_p = np.array(a_p, dtype=int)
            ia_p = np.array(ia_p, dtype=int)

            # Sort according to energy difference
            p_s = np.argsort(w_p)
            f_p = f_p[p_s]
            w_p = w_p[p_s]
            i_p = i_p[p_s]
            a_p = a_p[p_s]
            ia_p = ia_p[p_s]

            Np = len(f_p)
            P_p = []
            for p in range(Np):
                P = np.ravel_multi_index(ia_p[p], (Nn, Nn))
                P_p.append(P)
            P_p = np.array(P_p)

            dm_vp = np.empty((3, Np), dtype=float)

        for v in range(3):
            direction = np.zeros(3, dtype=float)
            direction[v] = 1.0
            cef = ConstantElectricField(Hartree / Bohr, direction)
            kick_hamiltonian = KickHamiltonian(paw.hamiltonian, paw.density,
                                               cef)
            dm_MM = paw.wfs.eigensolver.calculate_hamiltonian_matrix(
                kick_hamiltonian, paw.wfs, paw.wfs.kpt_u[u],
                add_kinetic=False, root=-1)

            if self.ksl.using_blacs:
                tmp_nM = self.ksl.mmdescriptor.zeros(dtype=C0_nM.dtype)
                pblas_simple_hemm(self.ksl.mmdescriptor,
                                  self.ksl.mmdescriptor,
                                  self.ksl.mmdescriptor,
                                  dm_MM, C0_nM.conj(), tmp_nM,
                                  side='R', uplo='L')
                dm_nn = self.ksl.mmdescriptor.zeros(dtype=C0_nM.dtype)
                pblas_simple_gemm(self.ksl.mmdescriptor,
                                  self.ksl.mmdescriptor,
                                  self.ksl.mmdescriptor,
                                  tmp_nM, C0_nM, dm_nn, transb='T')
            else:
                tri2full(dm_MM)
                dm_nn = np.dot(C0_nM.conj(), np.dot(dm_MM, C0_nM.T))

            dm_nn = collect_MM(self.ksl, dm_nn)
            if self.comm.rank == 0:
                dm_P = dm_nn.ravel()
                dm_p = dm_P[P_p]
                dm_vp[v] = dm_p

        if self.comm.rank == 0:
            self.w_p = w_p
            self.f_p = f_p
            self.ia_p = ia_p
            self.P_p = P_p
            self.dm_vp = dm_vp

        self.has_initialized = True

    def write(self, filename):
        from ase.io.trajectory import write_atoms

        self.log(f'{self.__class__.__name__}: Writing to {filename}')
        writer = Writer(filename, self.world, mode='w',
                        tag=self.__class__.ulmtag)
        writer.write(version=self.__class__.version)

        write_atoms(writer.child('atoms'), self.atoms)

        writer.write(ha=Hartree)
        write_uMM(self.kd, self.ksl, writer, 'S_uMM', self.S_uMM)
        write_uMM(self.kd, self.ksl, writer, 'C0_unM', self.C0_unM)
        write_uX(self.kd, self.ksl.block_comm, writer, 'eig_un', self.eig_un)
        write_uX(self.kd, self.ksl.block_comm, writer, 'occ_un', self.occ_un)

        if self.comm.rank == 0:
            for arg in self.readwrite_attrs:
                writer.write(arg, getattr(self, arg))

        writer.close()

    def read(self, filename):
        self.reader = Reader(filename)
        tag = self.reader.get_tag()
        if tag != self.__class__.ulmtag:
            raise RuntimeError('Unknown tag %s' % tag)
        self.version = self.reader.version

        # Do lazy reading in __getattr__ only if/when
        # the variables are required
        self.has_initialized = True

    def __getattr__(self, attr):
        if attr in ['S_uMM', 'C0_unM']:
            val = read_uMM(self.kpt_u, self.ksl, self.reader, attr)
            setattr(self, attr, val)
            return val
        if attr in ['eig_un', 'occ_un']:
            val = read_uX(self.kpt_u, self.reader, attr)
            setattr(self, attr, val)
            return val
        if attr in ['C0S_unM']:
            C0S_unM = []
            for u, kpt in enumerate(self.kpt_u):
                C0_nM = self.C0_unM[u]
                S_MM = self.S_uMM[u]
                if self.ksl.using_blacs:
                    C0S_nM = self.ksl.mmdescriptor.zeros(dtype=C0_nM.dtype)
                    pblas_simple_hemm(self.ksl.mmdescriptor,
                                      self.ksl.mmdescriptor,
                                      self.ksl.mmdescriptor,
                                      S_MM, C0_nM, C0S_nM,
                                      side='R', uplo='L')
                else:
                    C0S_nM = np.dot(C0_nM, S_MM)
                C0S_unM.append(C0S_nM)
            setattr(self, attr, C0S_unM)
            return C0S_unM
        if attr in ['weight_Mn']:
            assert self.world.size == 1
            C2_nM = np.absolute(self.C0_unM[0])**2
            val = C2_nM.T / np.sum(C2_nM, axis=1)
            setattr(self, attr, val)
            return val

        try:
            val = getattr(self.reader, attr)
            if attr == 'atoms':
                from ase.io.trajectory import read_atoms
                val = read_atoms(val)
            setattr(self, attr, val)
            return val
        except (KeyError, AttributeError):
            pass

        raise AttributeError('Attribute %s not defined in version %s' %
                             (repr(attr), repr(self.version)))

    def distribute(self, comm):
        self.comm = comm
        N = comm.size
        self.Np = len(self.P_p)
        self.Nq = int(np.ceil(self.Np / float(N)))
        self.NQ = self.Nq * N
        self.w_q = self.distribute_p(self.w_p)
        self.f_q = self.distribute_p(self.f_p)
        self.dm_vq = self.distribute_xp(self.dm_vp)

    def distribute_p(self, a_p, a_q=None, root=0):
        if a_q is None:
            a_q = np.zeros(self.Nq, dtype=a_p.dtype)
        if self.comm.rank == root:
            a_Q = np.append(a_p, np.zeros(self.NQ - self.Np, dtype=a_p.dtype))
        else:
            a_Q = None
        self.comm.scatter(a_Q, a_q, root)
        return a_q

    def collect_q(self, a_q, root=0):
        if self.comm.rank == root:
            a_Q = np.zeros(self.NQ, dtype=a_q.dtype)
        else:
            a_Q = None
        self.comm.gather(a_q, root, a_Q)
        if self.comm.rank == root:
            a_p = a_Q[:self.Np]
        else:
            a_p = None
        return a_p

    def distribute_xp(self, a_xp):
        Nx = a_xp.shape[0]
        a_xq = np.zeros((Nx, self.Nq), dtype=a_xp.dtype)
        for x in range(Nx):
            self.distribute_p(a_xp[x], a_xq[x])
        return a_xq

    def transform(self, rho_uMM, broadcast=False):
        assert len(rho_uMM) == 1, 'K-points not implemented'
        u = 0
        rho_MM = np.ascontiguousarray(rho_uMM[u])
        C0S_nM = self.C0S_unM[u].astype(rho_MM.dtype, copy=True)
        # KS decomposition
        if self.ksl.using_blacs:
            tmp_nM = self.ksl.mmdescriptor.zeros(dtype=rho_MM.dtype)
            pblas_simple_gemm(self.ksl.mmdescriptor,
                              self.ksl.mmdescriptor,
                              self.ksl.mmdescriptor,
                              C0S_nM, rho_MM, tmp_nM)
            rho_nn = self.ksl.mmdescriptor.zeros(dtype=rho_MM.dtype)
            pblas_simple_gemm(self.ksl.mmdescriptor,
                              self.ksl.mmdescriptor,
                              self.ksl.mmdescriptor,
                              tmp_nM, C0S_nM, rho_nn, transb='C')
        else:
            rho_nn = np.dot(np.dot(C0S_nM, rho_MM), C0S_nM.T.conj())

        rho_nn = collect_MM(self.ksl, rho_nn)
        if self.comm.rank == 0:
            rho_P = rho_nn.ravel()
            # Remove de-excitation terms
            rho_p = rho_P[self.P_p]
            if self.only_ia:
                rho_p *= 2
        else:
            rho_p = None

        if broadcast:
            if self.comm.rank != 0:
                rho_p = np.zeros_like(self.P_p, dtype=rho_MM.dtype)
            self.comm.broadcast(rho_p, 0)
        rho_up = [rho_p]
        return rho_up

    def ialims(self):
        i_p = self.ia_p[:, 0]
        a_p = self.ia_p[:, 1]
        imin = np.min(i_p)
        imax = np.max(i_p)
        amin = np.min(a_p)
        amax = np.max(a_p)
        return imin, imax, amin, amax

    def M_p_to_M_ia(self, M_p):
        return self.M_ia_from_M_p(M_p)

    def M_ia_from_M_p(self, M_p):
        imin, imax, amin, amax = self.ialims()
        M_ia = np.zeros((imax - imin + 1, amax - amin + 1), dtype=M_p.dtype)
        for M, (i, a) in zip(M_p, self.ia_p):
            M_ia[i - imin, a - amin] = M
        return M_ia

    def plot_matrix(self, M_p):
        import matplotlib.pyplot as plt
        M_ia = self.M_ia_from_M_p(M_p)
        plt.imshow(M_ia, interpolation='none')
        plt.xlabel('a')
        plt.ylabel('i')

    def get_dipole_moment_contributions(self, rho_up):
        assert len(rho_up) == 1, 'K-points not implemented'
        u = 0
        rho_p = rho_up[u]
        dmrho_vp = - self.dm_vp * rho_p
        return dmrho_vp

    def get_dipole_moment(self, rho_up):
        assert len(rho_up) == 1, 'K-points not implemented'
        u = 0
        rho_p = rho_up[u]
        dm_v = - np.dot(self.dm_vp, rho_p)
        return dm_v

    def get_density(self, wfs, rho_up, density='comp'):
        from gpaw.lcaotddft.densitymatrix import get_density

        if self.ksl.using_blacs:
            raise NotImplementedError('Scalapack is not supported')

        density_type = density
        # assert len(rho_up) == 1, 'K-points not implemented'
        u = 0
        #rho_p = rho_up[u]
        rho_p = rho_up
        C0_nM = self.C0_unM[u]

        rho_ia = self.M_ia_from_M_p(rho_p)
        imin, imax, amin, amax = self.ialims()
        C0_iM = C0_nM[imin:(imax + 1)]
        C0_aM = C0_nM[amin:(amax + 1)]

        rho_MM = np.dot(C0_iM.T, np.dot(rho_ia, C0_aM.conj()))
        rho_MM = 0.5 * (rho_MM + rho_MM.T)

        return get_density(rho_MM, wfs, self.density, density_type, u)

    def get_contributions_table(self, weight_p, minweight=0.01,
                                zero_fermilevel=True):
        assert weight_p.dtype == float
        u = 0  # TODO

        absweight_p = np.absolute(weight_p)
        tot_weight = weight_p.sum()
        propweight_p = weight_p / tot_weight * 100
        tot_propweight = propweight_p.sum()
        rest_weight = tot_weight
        rest_propweight = tot_propweight
        eig_n = self.eig_un[u].copy()
        if zero_fermilevel:
            eig_n -= self.fermilevel

        txt = ''
        txt += ('# %6s %4s(%8s)    %4s(%8s)  %12s %14s %8s\n' %
                ('p', 'i', 'eV', 'a', 'eV', 'Ediff (eV)', 'weight', '%'))
        p_s = np.argsort(absweight_p)[::-1]
        for s, p in enumerate(p_s):
            i, a = self.ia_p[p]
            if absweight_p[p] < minweight:
                break
            txt += ('  %6s %4d(%8.3f) -> %4d(%8.3f): %12.4f %14.4f %8.1f\n' %
                    (p, i, eig_n[i] * Hartree, a, eig_n[a] * Hartree,
                     self.w_p[p] * Hartree, weight_p[p], propweight_p[p]))
            rest_weight -= weight_p[p]
            rest_propweight -= propweight_p[p]
        txt += ('  %39s: %12s %+14.4f %8.1f\n' %
                ('rest', '', rest_weight, rest_propweight))
        txt += ('  %39s: %12s %+14.4f %8.1f\n' %
                ('total', '', tot_weight, tot_propweight))
        return txt

    def plot_TCM(self, weight_p, energy_o, energy_u, sigma,
                 zero_fermilevel=True, vmax='80%'):
        from gpaw.lcaotddft.tcm import TCMPlotter
        plotter = TCMPlotter(self, energy_o, energy_u, sigma, zero_fermilevel)
        ax_tcm = plotter.plot_TCM(weight_p, vmax)
        ax_occ_dos, ax_unocc_dos = plotter.plot_DOS()
        return ax_tcm, ax_occ_dos, ax_unocc_dos

    def get_TCM(self, weight_p, eig_n, energy_o, energy_u, sigma):
        flt_p = self.filter_by_x_ia(eig_n, energy_o, energy_u, 8 * sigma)
        weight_f = weight_p[flt_p]
        G_fo = gauss_ij(eig_n[self.ia_p[flt_p, 0]], energy_o, sigma)
        G_fu = gauss_ij(eig_n[self.ia_p[flt_p, 1]], energy_u, sigma)
        tcm_ou = np.dot(G_fo.T * weight_f, G_fu)
        return tcm_ou

    def get_DOS(self, eig_n, energy_o, energy_u, sigma):
        return self.get_weighted_DOS(1, eig_n, energy_o, energy_u, sigma)

    def get_weighted_DOS(self, weight_n, eig_n, energy_o, energy_u, sigma):
        if not isinstance(weight_n, np.ndarray):
            # Assume float
            weight_n = weight_n * np.ones_like(eig_n)
        G_on = gauss_ij(energy_o, eig_n, sigma)
        G_un = gauss_ij(energy_u, eig_n, sigma)
        dos_o = np.dot(G_on, weight_n)
        dos_u = np.dot(G_un, weight_n)
        return dos_o, dos_u

    def get_weight_n_by_l(self, l):
        if isinstance(l, int):
            weight_n = np.sum(self.weight_Mn[self.l_M == l], axis=0)
        else:
            weight_n = np.sum([self.get_weight_n_by_l(l_) for l_ in l],
                              axis=0)
        return weight_n

    def get_weight_n_by_a(self, a):
        if isinstance(a, int):
            weight_n = np.sum(self.weight_Mn[self.a_M == a], axis=0)
        else:
            weight_n = np.sum([self.get_weight_n_by_a(a_) for a_ in a],
                              axis=0)
        return weight_n

    def get_distribution_i(self, weight_p, energy_e, sigma,
                           zero_fermilevel=True):
        eig_n, fermilevel = self.get_eig_n(zero_fermilevel)
        flt_p = self.filter_by_x_i(eig_n, energy_e, 8 * sigma)
        weight_f = weight_p[flt_p]
        G_fe = gauss_ij(eig_n[self.ia_p[flt_p, 0]], energy_e, sigma)
        dist_e = np.dot(G_fe.T, weight_f)
        return dist_e

    def get_distribution_a(self, weight_p, energy_e, sigma,
                           zero_fermilevel=True):
        eig_n, fermilevel = self.get_eig_n(zero_fermilevel)
        flt_p = self.filter_by_x_a(eig_n, energy_e, 8 * sigma)
        weight_f = weight_p[flt_p]
        G_fe = gauss_ij(eig_n[self.ia_p[flt_p, 1]], energy_e, sigma)
        dist_e = np.dot(G_fe.T, weight_f)
        return dist_e

    def get_distribution_ia(self, weight_p, energy_o, energy_u, sigma,
                            zero_fermilevel=True):
        """
        Filter both i and a spaces as in TCM.

        """
        eig_n, fermilevel = self.get_eig_n(zero_fermilevel)
        flt_p = self.filter_by_x_ia(eig_n, energy_o, energy_u, 8 * sigma)
        weight_f = weight_p[flt_p]
        G_fo = gauss_ij(eig_n[self.ia_p[flt_p, 0]], energy_o, sigma)
        dist_o = np.dot(G_fo.T, weight_f)
        G_fu = gauss_ij(eig_n[self.ia_p[flt_p, 1]], energy_u, sigma)
        dist_u = np.dot(G_fu.T, weight_f)
        return dist_o, dist_u

    def get_distribution(self, weight_p, energy_e, sigma):
        w_p = self.w_p * Hartree
        flt_p = self.filter_by_x_p(w_p, energy_e, 8 * sigma)
        weight_f = weight_p[flt_p]
        G_fe = gauss_ij(w_p[flt_p], energy_e, sigma)
        dist_e = np.dot(G_fe.T, weight_f)
        return dist_e

    def get_eig_n(self, zero_fermilevel=True):
        u = 0  # TODO
        eig_n = self.eig_un[u].copy()
        if zero_fermilevel:
            eig_n -= self.fermilevel
            fermilevel = 0.0
        else:
            fermilevel = self.fermilevel
        eig_n *= Hartree
        fermilevel *= Hartree
        return eig_n, fermilevel

    def filter_by_x_p(self, x_p, energy_e, buf):
        flt_p = np.logical_and((energy_e[0] - buf) <= x_p,
                               x_p <= (energy_e[-1] + buf))
        return flt_p

    def filter_by_x_i(self, x_n, energy_e, buf):
        return self.filter_by_x_p(x_n[self.ia_p[:, 0]], energy_e, buf)

    def filter_by_x_a(self, x_n, energy_e, buf):
        return self.filter_by_x_p(x_n[self.ia_p[:, 1]], energy_e, buf)

    def filter_by_x_ia(self, x_n, energy_o, energy_u, buf):
        flti_p = self.filter_by_x_i(x_n, energy_o, buf)
        flta_p = self.filter_by_x_a(x_n, energy_u, buf)
        flt_p = np.logical_and(flti_p, flta_p)
        return flt_p

    def __del__(self):
        if self.reader is not None:
            self.reader.close()
