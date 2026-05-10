import numpy as np
from guesthost.analysis.functions_modules import (
        proj2plane, self_correlation_P1, self_correlation_P2,
        distance_of_point_on_plane
)
from pymatgen.util.coord_cython import pbc_shortest_vectors
from pymatgen.core import Lattice as PmgLattice

class Motif:
    """Class representing one unit of a repeating collection of fragments"""
    
    def __init__(self, frg_list, unitcell_data=None, br_cage_coords=None, br_cage_indices=None):
        self.fragments = frg_list
        self.nfrg = len(frg_list)
        self.unitcell_data = unitcell_data
        self.br_cage_coords = br_cage_coords
        self.br_cage_indices = br_cage_indices

    @property
    def ma(self):
        """Methylammonium fragment for hybrid perovskite motifs."""
        return self.fragments[0]

    @property
    def host(self):
        """Host fragment for hybrid perovskite motifs."""
        return self.fragments[1]

class Lattice:
    """Class representing a lattice composed of several motifs with its associated cell"""

    def __init__(self, motif_list, cell, nx, ny, nz, pbc=True):
        motif_grid = np.array(motif_list).reshape(nx, ny, nz)

        self.motif_grid = motif_grid
        self.cell = cell
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.ncells = (nx, ny, nz)
        self.pbc = pbc

    def compute_for_all(self, fn, *args, dtype=float, **kwargs):
        all_vals = np.zeros((self.nx, self.ny, self.nz), dtype=dtype)

        for ix, iy, iz in np.ndindex(self.nx, self.ny, self.nz):
            val = fn(self, (ix, iy, iz), *args, **kwargs)

            all_vals[ix, iy, iz] = val

        return all_vals


class HPLattice(Lattice):
    """Class representing a hybrid perovskite lattice"""

    def __init__(self, motif_list, cell, nx, ny, nz):
        super().__init__(motif_list, cell, nx, ny, nz)

    def ma_torsion(self, ind, frg_ref, htyp="N", atan=False):
        frg = self.get_motif(ind).ma
        if htyp == "N":
            v1 = frg_ref.HN
        elif htyp == "C":
            v1 = frg_ref.HC

        if atan:
            val = frg.caltorsion_atan(v1, htyp=htyp)
        else:
            val = frg.caltorsion(v1, htyp=htyp)

        return val

    def all_ma_torsions(self, lattice_ref, htyp="N", atan=False):
        torsions = np.zeros((self.nx, self.ny, self.nz))
        for ix, iy, iz in np.ndindex(self.nx, self.ny, self.nz):
            frg_ref = lattice_ref.get_motif((ix, iy, iz)).ma
            torsion = self.ma_torsion((ix, iy, iz), frg_ref, htyp=htyp, atan=atan)

            torsions[ix,iy,iz] = torsion

        return torsions

    def ma_torsion_removetilts(self, ind, frg_ref, htyp="N", match_rows=False, return_tilt=False):
        frg = self.get_motif(ind).ma
        if htyp == "N":
            v1 = frg_ref.HN
        elif htyp == "C":
            v1 = frg_ref.HC

        val = frg.caltorsion_removetilts(
            v1, htyp=htyp, match_rows=match_rows, 
            return_tilt=return_tilt
        )

        return val

    def all_ma_torsions_removetilts(self, lattice_ref, htyp="N", match_rows=False):
        torsions = np.zeros((self.nx, self.ny, self.nz))

        for ix, iy, iz in np.ndindex(self.nx, self.ny, self.nz):
            frg_ref = lattice_ref.get_motif((ix, iy, iz)).ma
            torsion = self.ma_torsion_removetilts(
                    (ix, iy, iz), frg_ref, htyp=htyp, 
                    match_rows=match_rows, return_tilt=False
                    )

            torsions[ix,iy,iz] = torsion

        return torsions

    def all_ma_torsions_withtilts(self, lattice_ref, htyp="N", match_rows=False):
        torsions = np.zeros((self.nx, self.ny, self.nz))
        tilts = np.zeros((self.nx, self.ny, self.nz))

        for ix, iy, iz in np.ndindex(self.nx, self.ny, self.nz):
            frg_ref = lattice_ref.get_motif((ix, iy, iz)).ma
            torsion, tilt = self.ma_torsion_removetilts(
                    (ix, iy, iz), frg_ref, htyp=htyp, 
                    match_rows=match_rows, return_tilt=True
                    )

            torsions[ix,iy,iz] = torsion
            tilts[ix,iy,iz] = tilt

        return torsions, tilts

    def ucell_localcellvecs(self, ind):
        """Compute local cell vectors from the host fragments around a motif."""
        return [np.array(v) for v in np.array(self.pb_environment_cell(ind)).T]

    def _localcellparameters_from_vecs(self, vecs):
        av, bv, cv = vecs
        a, b, c = [np.linalg.norm(v) for v in [av, bv, cv]]

        alpha = np.degrees(np.arccos(np.clip(np.dot(bv, cv) / (b * c), -1.0, 1.0)))
        beta = np.degrees(np.arccos(np.clip(np.dot(av, cv) / (a * c), -1.0, 1.0)))
        gamma = np.degrees(np.arccos(np.clip(np.dot(av, bv) / (a * b), -1.0, 1.0)))

        cellmat = np.column_stack([av, bv, cv])
        vol = np.linalg.det(cellmat)

        return {
            'a': a, 'b': b, 'c': c,
            'alpha': alpha, 'beta': beta, 'gamma': gamma,
            'α': alpha, 'β': beta, 'γ': gamma,
            'vol': vol,
        }

    def ucell_localcellparameters(self, ind):
        """Compute local cell parameters from the host fragments around a motif."""
        return self._localcellparameters_from_vecs(self.ucell_localcellvecs(ind))

    def ucell_localcelllengths(self, ind):
        """Compute local cell vector lengths."""
        return np.array([np.linalg.norm(v) for v in self.ucell_localcellvecs(ind)])

    def ucell_localcelllengths_only(self, ind):
        """Compute local cell vector lengths as a list."""
        return self.ucell_localcelllengths(ind).tolist()

    def ucell_localcellangles(self, ind):
        """Compute local cell angles."""
        params = self.ucell_localcellparameters(ind)
        return [params["alpha"], params["beta"], params["gamma"]]

    def ucell_localcell_lengths_angles(self, ind):
        """Compute local cell lengths and angles."""
        return self.ucell_localcelllengths_only(ind), self.ucell_localcellangles(ind)

    def _ortho_axis_indices(self, ind, dir, small=False):
        diagonal_shifts = (
            ((0, 1, -1), (0, 1, 1)),
            ((-1, 0, 1), (1, 0, 1)),
            ((1, 1, 0), (-1, 1, 0)),
        )
        base_shift = [0, 0, 0]
        base_shift[dir] = 1 if small else 2
        shifts = [tuple(base_shift), *diagonal_shifts[dir]]
        return [self.shift_index_vector(ind, shift) for shift in shifts]

    def ucell_localcellvecs_ortho(self, ind, dir, small=False):
        """Compute orthorhombic-like local cell vectors from Pb motifs."""
        origin = self.get_motif(ind).host.B
        vecs = []
        for shifted_ind in self._ortho_axis_indices(ind, dir, small=small):
            pos = self.get_motif(shifted_ind).host.B
            if self.pbc:
                vec = pbc_vectors(self.cell, origin, pos)[0]
            else:
                vec = pos - origin
            vecs.append(vec)
        if small:
            vecs[0] = 2 * vecs[0]
        return vecs

    def ucell_localcelllengths_ortho(self, ind, dir, small=False):
        """Compute orthorhombic-like local cell vector lengths."""
        return [np.linalg.norm(v) for v in self.ucell_localcellvecs_ortho(ind, dir=dir, small=small)]

    def ucell_localcellangles_ortho(self, ind, dir, small=False):
        """Compute orthorhombic-like local cell angles."""
        av, bv, cv = self.ucell_localcellvecs_ortho(ind, dir=dir, small=small)
        return [
            self._angle_between(bv, cv),
            self._angle_between(av, cv),
            self._angle_between(av, bv),
        ]

    def ucell_localcell_lengths_angles_ortho(self, ind, dir, small=False):
        """Compute orthorhombic-like local cell lengths and angles."""
        return (
            self.ucell_localcelllengths_ortho(ind, dir=dir, small=small),
            self.ucell_localcellangles_ortho(ind, dir=dir, small=small),
        )

    def _angle_between(self, v1, v2):
        n1 = np.linalg.norm(v1)
        n2 = np.linalg.norm(v2)
        if n1 == 0 or n2 == 0:
            return 0.0
        return np.degrees(np.arccos(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0)))

    def reciprocalvecs(self, h):
        """Compute reciprocal vectors for cell matrix h"""
        v = np.cross(h[1], h[2])
        vol = np.dot(h[0], v)
        g = np.array([np.cross(h[j], h[k]) / vol for i, j, k in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]])
        return g

    def _theta_phi_from_vector(self, r, lcell, dir):
        """Compute theta and phi for a vector relative to local cell row vectors."""
        lst = [(0, 1, 2), (1, 2, 0), (2, 0, 1)]
        i, j, k = lst[dir]

        lcell = np.array(lcell)
        theta = np.degrees(np.arccos(np.clip(
            np.dot(r, lcell[i]) / (np.linalg.norm(r) * np.linalg.norm(lcell[i])),
            -1.0,
            1.0,
        )))

        g = self.reciprocalvecs(lcell)
        vprp = np.dot(r, g[i]) * lcell[i]
        vpl = r - vprp
        vpl_norm = np.linalg.norm(vpl)
        if vpl_norm < 1e-12:
            return theta, 90.0

        angb = np.dot(vpl, lcell[j]) / (vpl_norm * np.linalg.norm(lcell[j]))
        angc = np.dot(vpl, g[k])
        phase = 1.0 if angc >= 0 else -1.0
        phi = phase * np.degrees(np.arccos(np.clip(angb, -1, 1)))
        if np.isnan(phi):
            phi = 90.0

        return theta, phi

    def ucell_theta_phi(self, ind, dir):
        """Compute theta and phi for the MA fragment in a motif."""
        ma = self.get_motif(ind).ma
        return self._theta_phi_from_vector(ma.N - ma.C, self.pb_environment_cell(ind), dir)

    def ucell_theta_phi_dict(self, ind, dir):
        """Compute theta and phi as a mapping."""
        theta, phi = self.ucell_theta_phi(ind, dir)
        return {"theta": theta, "phi": phi}

    def ucell_theta_phi_withcutoff(self, ind, dir, ratio_cut=0.05):
        """Compute theta/phi and set phi to zero when the in-plane component is small."""
        ma = self.get_motif(ind).ma
        r = ma.N - ma.C
        lcell = np.array(self.pb_environment_cell(ind))
        theta, phi = self._theta_phi_from_vector(r, lcell, dir)

        i, j, k = [(0, 1, 2), (1, 2, 0), (2, 0, 1)][dir]
        g = self.reciprocalvecs(lcell)
        vprp = np.dot(r, g[i]) * lcell[i]
        vpl = r - vprp
        cratio = np.linalg.norm(vpl) / np.linalg.norm(r)
        if cratio <= ratio_cut:
            phi = 0.0
        else:
            vpl_norm = np.linalg.norm(vpl)
            angb = np.dot(vpl, lcell[j]) / (vpl_norm * np.linalg.norm(lcell[j]))
            angc = np.dot(vpl, g[k])
            phase = -1.0 if angc < 0 else 1.0
            phi = phase * np.degrees(np.arccos(np.clip(angb, -1.0, 1.0)))
        return {"theta": theta, "phi": phi}

    @staticmethod
    def ucell_xyz(theta_deg, phi_deg, theta_dir=2):
        """Convert theta/phi angles to xyz components."""
        order_map = [[2, 0, 1], [1, 2, 0], [0, 1, 2]]
        theta = np.deg2rad(theta_deg)
        phi = np.deg2rad(phi_deg)
        vec = [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ]
        return [vec[i] for i in order_map[theta_dir]]

    def all_ucell_theta_phi(self, dir):
        """Compute MA theta/phi for every motif in the lattice."""
        theta_vals = np.zeros((self.nx, self.ny, self.nz))
        phi_vals = np.zeros((self.nx, self.ny, self.nz))

        for ind in np.ndindex(self.nx, self.ny, self.nz):
            theta, phi = self.ucell_theta_phi(ind, dir)
            theta_vals[ind] = theta
            phi_vals[ind] = phi

        return theta_vals, phi_vals

    def ucell_br_distance_from_plane(self, ind, dir_pl, dir_br):
        """Compute Br distance from a Pb plane using host fragments."""
        return self.br_distance_from_plane(ind, dir_pl, dir_br)

    def ucell_br_br_vectors(self, ind, dir=0):
        """Compute Br(-axis) to Br(+axis) vectors for a motif."""
        axes_order = [[0, 1, 2], [2, 0, 1], [1, 2, 0]][dir]
        host = self.get_motif(ind).host
        vecs = []
        for ax in axes_order:
            start = self.get_motif(self.shift_index(ind, ax, -1)).host.X[ax]
            end = host.X[ax]
            if self.pbc:
                vec = pbc_vectors(self.cell, start, end)[0]
            else:
                vec = end - start
            vecs.append(vec)
        return vecs

    def ucell_br_br_distances(self, ind, dir=0):
        """Compute Br(-axis) to Br(+axis) distances."""
        return [np.linalg.norm(v) for v in self.ucell_br_br_vectors(ind, dir=dir)]

    def ucell_tetragonal_lengthening_distortion_br(self, ind, dir=0):
        """Compute Br cage tetragonal lengthening distortion."""
        br_dists = self.ucell_br_br_distances(ind, dir=dir)
        return (2 * br_dists[0] / (br_dists[1] + br_dists[2])) - 1

    def ucell_scissoring_distortion_br(self, ind, dir=0):
        """Compute Br cage scissoring distortion."""
        br_vecs = self.ucell_br_br_vectors(ind, dir=dir)
        v2 = br_vecs[1] / np.linalg.norm(br_vecs[1])
        v3 = br_vecs[2] / np.linalg.norm(br_vecs[2])
        return float(np.dot(v2, v3))

    def ucell_scissoring_distortion_angle_br(self, ind, dir=0):
        """Compute Br cage scissoring distortion angle."""
        omega = self.ucell_scissoring_distortion_br(ind, dir=dir)
        return 90.0 - np.degrees(np.arccos(np.clip(omega, -1.0, 1.0)))

    def ucell_pb_br_pb_angle(self, ind, dir):
        """Compute Pb-Br-Pb angle using host fragments."""
        return self.pb_br_pb_angle(ind, dir)

    def ucell_pb_br_pb_angle_hostrelative(self, ind, dir, dir_coup):
        """Compute host-relative Pb-Br-Pb angle using host fragments."""
        return self.pb_br_pb_angle_hostrelative(ind, dir, dir_coup)

    def ucell_pb_br_pb_angle_alongaxis(self, ind, dir, dir_coup):
        """Compute Pb-Br-Pb angle after projection along the coupled host axis."""
        o_host = self.get_motif(ind).host
        ax_hosts = [self.get_motif(self.shift_index(ind, ax)).host for ax in range(3)]
        perp_dirs = [(1, 2), (2, 0), (0, 1)][dir_coup]

        if self.pbc:
            v1 = pbc_vectors(self.cell, o_host.B, ax_hosts[perp_dirs[0]].B)[0]
            v2 = pbc_vectors(self.cell, o_host.B, ax_hosts[perp_dirs[1]].B)[0]
            vAB = pbc_vectors(self.cell, o_host.X[dir], o_host.B)[0]
            vCB = pbc_vectors(self.cell, o_host.X[dir], ax_hosts[dir].B)[0]
        else:
            v1 = ax_hosts[perp_dirs[0]].B - o_host.B
            v2 = ax_hosts[perp_dirs[1]].B - o_host.B
            vAB = o_host.B - o_host.X[dir]
            vCB = ax_hosts[dir].B - o_host.X[dir]

        axis = np.cross(v1, v2)
        n = axis / np.linalg.norm(axis)
        vABp = vAB - np.dot(vAB, n) * n
        vCBp = vCB - np.dot(vCB, n) * n
        normABp = np.linalg.norm(vABp)
        normCBp = np.linalg.norm(vCBp)
        if normABp < 1e-12 or normCBp < 1e-12:
            return 0.0
        costheta = np.clip(np.dot(vABp, vCBp) / (normABp * normCBp), -1.0, 1.0)
        return np.degrees(np.arccos(costheta))

    def ucell_pb_br_pb_angle_alongaxis_hostrelative(self, ind, dir, dir_coup):
        """Compute host-relative projected Pb-Br-Pb angle."""
        ang = self.ucell_pb_br_pb_angle_alongaxis(ind, dir, dir_coup)
        pl_dist = self.ucell_br_distance_from_plane(ind, dir_pl=[dir, dir_coup], dir_br=dir)
        return 360.0 - ang if pl_dist > 0 else ang

    def ucell_η_lat(self, ind, dir_coup, ref=180.0, N=1.0):
        """Compute lattice order parameter using host fragments."""
        br_dirs = [(1, 2), (2, 0), (0, 1)][dir_coup]
        θ_v = [self.ucell_pb_br_pb_angle_hostrelative(ind, d, dir_coup) for d in br_dirs]
        V = (np.array(θ_v) - ref) / N
        return tuple(V.tolist())

    ucell_eta_lat = ucell_η_lat

    def all_ucell_η_lat(self, dir_coup, ref=180.0, N=1.0):
        vals = np.zeros((self.nx, self.ny, self.nz, 2))
        for ind in np.ndindex(self.nx, self.ny, self.nz):
            vals[ind] = self.ucell_η_lat(ind, dir_coup, ref=ref, N=N)
        return vals

    all_ucell_eta_lat = all_ucell_η_lat

    def ucell_η_lat_alongaxis(self, ind, dir_coup, ref=180.0, N=1.0):
        """Compute lattice order parameter from projected Pb-Br-Pb angles."""
        br_dirs = [(1, 2), (2, 0), (0, 1)][dir_coup]
        θ_v = [
            self.ucell_pb_br_pb_angle_alongaxis_hostrelative(ind, d, dir_coup)
            for d in br_dirs
        ]
        V = (np.array(θ_v) - ref) / N
        return tuple(V.tolist())

    ucell_eta_lat_alongaxis = ucell_η_lat_alongaxis

    def U_r(self, R, k_vec, MA=False):
        """Rotation matrix for order parameters"""
        k_dot_R = np.dot(k_vec, R)
        if MA:
            k_dot_R = (np.pi / 4) * (1.0 - np.cos(k_dot_R))
        c = np.cos(k_dot_R)
        s = np.sin(k_dot_R)
        U = np.array([[c, s], [-s, c]])
        return U

    def ucell_ω(self, η, R, k_vec=None, orig_cell=(0, 0), MA=False):
        """Transform η to ω using rotation matrix"""
        if k_vec is None:
            k_vec = np.array([0.5, 0.5]) * 2 * np.pi
        U = self.U_r(np.array(R) + np.array(orig_cell), k_vec, MA)
        ω = U @ η
        return np.array(ω)

    def _change_ϕ_domain(self, ϕ, orig_domain=(-180, 180), final_domain=(-45, 135), eps=None):
        """Change phi domain"""
        oi, of_ = orig_domain
        fi, ff = final_domain
        period = ff - fi
        if eps is not None and abs(ϕ) < eps + fi:
            ϕ = fi
        elif eps is not None and abs(ϕ) > ff - eps:
            ϕ = ff
        else:
            if oi <= ϕ <= fi:
                ϕ += period
            elif ff <= ϕ <= of_:
                ϕ -= period
        return ϕ

    def ucell_η_MA_static(self, θ, ϕ, norm=False, change_phi_domain=True, phi_domain=(-45, 135), eps=None):
        """Compute MA η from θ, ϕ"""
        if change_phi_domain:
            ϕ = self._change_ϕ_domain(ϕ, final_domain=phi_domain, eps=eps)
        vb = np.sin(np.radians(θ)) * np.cos(np.radians(ϕ))
        vc = np.sin(np.radians(θ)) * np.sin(np.radians(ϕ))
        V = np.array([vb, vc])
        if norm:
            n = np.linalg.norm(V)
            if n > 0:
                V /= n
        return V.tolist()

    ucell_eta_MA_from_angles = ucell_η_MA_static

    def ucell_η_MA(self, ind, dir, **kwargs):
        """Compute MA η from the MA fragment in a motif."""
        θ, ϕ = self.ucell_theta_phi(ind, dir)
        return self.ucell_η_MA_static(θ, ϕ, **kwargs)

    ucell_eta_MA = ucell_η_MA

    def all_ucell_η_MA(self, dir, **kwargs):
        vals = np.zeros((self.nx, self.ny, self.nz, 2))
        for ind in np.ndindex(self.nx, self.ny, self.nz):
            vals[ind] = self.ucell_η_MA(ind, dir, **kwargs)
        return vals

    all_ucell_eta_MA = all_ucell_η_MA

    def _change_θ_domain(self, θ, orig_domain=(0, 180), final_domain=(-45, 135)):
        """Change theta domain."""
        oi, of_ = orig_domain
        fi, ff = final_domain
        period = ff - fi
        if oi <= θ <= fi:
            θ += period
        elif ff <= θ <= of_:
            θ -= period
        return θ

    def ucell_η_3d_MA_static(
        self,
        θ,
        ϕ,
        norm=False,
        change_phi_domain=True,
        phi_domain=(-45, 135),
        eps=None,
        change_theta_domain=True,
        theta_domain=(-45, 135),
    ):
        """Compute 3D MA η from θ, ϕ."""
        if change_phi_domain:
            ϕ = self._change_ϕ_domain(ϕ, final_domain=phi_domain, eps=eps)
        if change_theta_domain:
            θ = self._change_θ_domain(θ, final_domain=theta_domain)
        va = np.sin(np.radians(θ)) * np.cos(np.radians(ϕ))
        vb = np.sin(np.radians(θ)) * np.sin(np.radians(ϕ))
        vc = np.cos(np.radians(θ))
        V = np.array([va, vb, vc])
        if norm:
            n = np.linalg.norm(V)
            if n > 0:
                V /= n
        return V.tolist()

    ucell_eta_3d_MA_from_angles = ucell_η_3d_MA_static

    def ucell_η_3d_MA(self, ind, dir=2, **kwargs):
        """Compute 3D MA η from the MA fragment in a motif."""
        θ, ϕ = self.ucell_theta_phi(ind, dir)
        return self.ucell_η_3d_MA_static(θ, ϕ, **kwargs)

    ucell_eta_3d_MA = ucell_η_3d_MA

    def ucell_MA_dihedrals(self, ind):
        """Compute all C-H-C-N-H-N MA dihedral angles with PBC."""
        ma = self.get_motif(ind).ma
        res = []
        for hn in ma.HN:
            for hc in ma.HC:
                res.append(_dihedral_pbc(self.cell, hc, ma.C, ma.N, hn))
        return tuple(res)

    def ucell_ma_displacement_vecs_br(self, ind):
        """Compute MA center displacement relative to the Br cage center."""
        motif = self.get_motif(ind)
        if motif.br_cage_coords is None:
            raise ValueError("Br cage coordinates are required for ucell_ma_displacement_vecs_br.")

        ma = motif.ma
        br_coords = np.array(motif.br_cage_coords)
        base = br_coords[0]
        if self.pbc:
            rel_br = pbc_vectors(self.cell, base, br_coords)
            rel_c = pbc_vectors(self.cell, base, ma.C)[0]
            rel_n = pbc_vectors(self.cell, base, ma.N)[0]
        else:
            rel_br = br_coords - base
            rel_c = ma.C - base
            rel_n = ma.N - base

        br_com = rel_br.mean(axis=0)
        ma_com = np.array([rel_c, rel_n]).mean(axis=0)
        if self.pbc:
            return pbc_vectors(self.cell, br_com, ma_com)[0]
        return ma_com - br_com

    def ucell_pb_axis(self, ind, python_index=False):
        """Return Pb axis indices stored by the motif model."""
        motif = self.get_motif(ind)
        origin = int(motif.host.indices[0])
        axes = [int(self.get_motif(self.shift_index(ind, ax)).host.indices[0]) for ax in range(3)]
        return (origin, axes)

    def ucell_udata(self, ind):
        """Return unit-cell index metadata reconstructed from motif fragments."""
        motif = self.get_motif(ind)
        if motif.unitcell_data is not None:
            return motif.unitcell_data
        ma = motif.ma
        host = motif.host
        br_axis_neg = [
            int(self.get_motif(self.shift_index(ind, ax, -1)).host.indices[1 + ax])
            for ax in range(3)
        ]
        return {
            "c_ind": int(ma.indices[0]),
            "n_ind": int(ma.indices[4]),
            "h_inds_c": [int(i) for i in ma.indices[1:4]],
            "h_inds_n": [int(i) for i in ma.indices[5:8]],
            "pb_inds": [int(self.get_motif(self.shift_index(ind, dx, sx)).host.indices[0])
                        for dx in range(3) for sx in (0, 1)],
            "br_inds": [] if motif.br_cage_indices is None else [int(i) for i in motif.br_cage_indices],
            "pb_axis": self.ucell_pb_axis(ind),
            "br_axis": [int(i) for i in host.indices[1:4]],
            "br_axis_neg": br_axis_neg,
        }

    def ucell_CN_distance(self, ind):
        """Compute the C-N distance in an MA fragment."""
        ma = self.get_motif(ind).ma
        if self.pbc:
            vec = pbc_vectors(self.cell, ma.C, ma.N)[0]
        else:
            vec = ma.N - ma.C
        return np.linalg.norm(vec)

    def ucell_CH_distances(self, ind):
        """Compute C-H distances for C-side hydrogens."""
        ma = self.get_motif(ind).ma
        if self.pbc:
            return tuple(np.linalg.norm(pbc_vectors(self.cell, ma.C, h)[0]) for h in ma.HC)
        return tuple(np.linalg.norm(h - ma.C) for h in ma.HC)

    def ucell_NH_distances(self, ind):
        """Compute N-H distances for N-side hydrogens."""
        ma = self.get_motif(ind).ma
        if self.pbc:
            return tuple(np.linalg.norm(pbc_vectors(self.cell, ma.N, h)[0]) for h in ma.HN)
        return tuple(np.linalg.norm(h - ma.N) for h in ma.HN)

    def ucell_all_octahedra_indices(self, ind):
        """Return the three Pb-Br-Pb octahedral axis index triplets."""
        origin, end_pbs = self.ucell_pb_axis(ind)
        mid_brs = [int(i) for i in self.get_motif(ind).host.indices[1:4]]
        return tuple((origin, mid_brs[i], end_pbs[i]) for i in range(3))

    def _pb_index_scaled_positions(self):
        cell = np.array(self.cell, dtype=float)
        invcell = np.linalg.inv(cell)
        entries = []
        for ind in np.ndindex(self.nx, self.ny, self.nz):
            host = self.get_motif(ind).host
            scaled = np.asarray(host.B, dtype=float) @ invcell
            entries.append((int(host.indices[0]), scaled, np.asarray(host.B, dtype=float)))
        return entries

    def get_global_index_data(self):
        """Return Pb index data used for orthorhombic global cell vectors."""
        entries = self._pb_index_scaled_positions()
        cell = np.array(self.cell, dtype=float)

        def closest_atom(target_frac):
            best = None
            best_d = np.inf
            target_frac = np.asarray(target_frac, dtype=float)
            for idx, scaled, _pos in entries:
                diff = scaled - target_frac
                diff -= np.round(diff)
                disp = diff @ cell
                dist = np.linalg.norm(disp)
                if dist < best_d:
                    best = idx
                    best_d = dist
            return best

        axes_OR = []
        for dir_coup in range(3):
            perp = [(1, 2), (2, 0), (0, 1)][dir_coup]
            fr_origin = np.zeros(3)
            fr_origin[perp[0]] = 0.5
            pb_origin = closest_atom(fr_origin)

            fr_endpt = np.zeros(3)
            fr_endpt[perp[1]] = 0.5
            pb_endpt = closest_atom(fr_endpt)

            axes_OR.append((pb_origin, pb_endpt))

        return {"supercell_size": self.nx, "axes_OR": axes_OR}

    def _scaled_by_pb_index(self):
        return {idx: scaled for idx, scaled, _pos in self._pb_index_scaled_positions()}

    def global_cell_vectors_OR(self, dir=0, global_data=None):
        """Compute global orthorhombic cell vectors from Pb positions."""
        if global_data is None:
            global_data = self.get_global_index_data()
        perp = [(1, 2), (2, 0), (0, 1)][dir]
        cell = np.array(self.cell, dtype=float)
        scaled_by_index = self._scaled_by_pb_index()
        a_OR = cell[dir]

        pb_origin, pb_endpt = global_data["axes_OR"][dir]
        fr_origin = scaled_by_index[pb_origin].copy()
        if fr_origin[dir] > 0.5:
            fr_origin[dir] -= 1.0
        if fr_origin[perp[1]] > 0.5:
            fr_origin[perp[1]] -= 1.0
        pos_origin = fr_origin @ cell

        fr_endpt = scaled_by_index[pb_endpt].copy()
        if fr_endpt[dir] > 0.5:
            fr_endpt[dir] -= 1.0
        if fr_endpt[perp[0]] > 0.5:
            fr_endpt[perp[0]] -= 1.0
        pos_endpt_c = fr_endpt @ cell

        fr_endpt_b = fr_endpt.copy()
        fr_endpt_b[perp[0]] += 1.0
        pos_endpt_b = fr_endpt_b @ cell

        b_OR = pos_endpt_b - pos_origin
        c_OR = pos_endpt_c - pos_origin
        return [np.array(a_OR), np.array(b_OR), np.array(c_OR)]

    def global_cell_vectors_OR_fromcell(self, dir=0):
        """Compute ideal orthorhombic cell vectors from the simulation cell."""
        dirb, dirc = [(1, 2), (2, 0), (0, 1)][dir]
        cell = np.array(self.cell, dtype=float)
        a_OR = 2.0 * cell[dir]
        b_OR = cell[dirb] - cell[dirc]
        c_OR = cell[dirb] + cell[dirc]
        return [a_OR, b_OR, c_OR]

    def global_cell_lengths_angles_OR(self, dir=0, global_data=None, fromcell=False):
        """Return global orthorhombic cell lengths and angles."""
        if fromcell:
            vecs = self.global_cell_vectors_OR_fromcell(dir=dir)
        else:
            vecs = self.global_cell_vectors_OR(dir=dir, global_data=global_data)
        params = self._localcellparameters_from_vecs(vecs)
        return (
            (params["a"], params["b"], params["c"]),
            (params["alpha"], params["beta"], params["gamma"]),
        )

    def global_cell_lengths_OR(self, dir=0, global_data=None, fromcell=False):
        """Return global orthorhombic cell lengths."""
        lengths, _angles = self.global_cell_lengths_angles_OR(
            dir=dir, global_data=global_data, fromcell=fromcell
        )
        return list(lengths)

    def global_cell_angles_OR(self, dir=0, global_data=None, fromcell=False):
        """Return global orthorhombic cell angles."""
        _lengths, angles = self.global_cell_lengths_angles_OR(
            dir=dir, global_data=global_data, fromcell=fromcell
        )
        return list(angles)

    def ma_orientation(self, ind, ax):
        "MA orientations for unitcell specified by `ind` with respect to axis, `ax`"""

        cn_vec = self.get_motif(ind).ma.CN
        lcell = self.pb_environment_cell(ind)
        vpl, t, p = proj2plane(ax, cn_vec, lcell)

        return t, p

    def all_ma_orientations(self, ax):
        theta_vals = np.zeros((self.nx, self.ny, self.nz))
        phi_vals = np.zeros((self.nx, self.ny, self.nz))

        for ix, iy, iz in np.ndindex(self.nx, self.ny, self.nz):
            t, p = self.ma_orientation((ix, iy, iz), ax)

            theta_vals[ix, iy, iz] = t
            phi_vals[ix, iy, iz] = p

        return theta_vals, phi_vals

    def shift_index(self, ind, ax, n=1):
        ind = list(ind)
        ind[ax] = (ind[ax]+n) % self.ncells[ax]

        return ind

    def shift_index_vector(self, ind, delta):
        return tuple((ind[ax] + delta[ax]) % self.ncells[ax] for ax in range(3))

    def get_motif(self, ind):
        ix, iy, iz = ind

        return self.motif_grid[ix, iy, iz]

    def pb_environment_coords(self, ind):
        ax_positions = np.array([
                self.get_motif(self.shift_index(ind, 0)).host.B,
                self.get_motif(self.shift_index(ind, 1)).host.B,
                self.get_motif(self.shift_index(ind, 2)).host.B
            ])

        return ax_positions

    def pb_environment_cell(self, ind):
        ix, iy, iz = ind
        o_host_frg = self.get_motif(ind).host
        o_pos = o_host_frg.B
        
        ax_positions = self.pb_environment_coords(ind)

        if self.pbc:
            lcell = pbc_vectors(self.cell, o_pos, ax_positions).T
        else:
            lcell = (ax_positions - o_pos).T

        return lcell

    def pb_br_pb_angle(self, ind, ax):
        ix, iy, iz = ind
        hostfrg_orig = self.get_motif(ind).host
        hostfrg_ax = self.get_motif(self.shift_index(ind, ax)).host

        atoms_orig = hostfrg_orig.atoms
        atoms_ax = hostfrg_ax.atoms
        atoms = atoms_orig + atoms_ax

        x_ind = ax+1
        ang = atoms.get_angle(0, x_ind, 4, mic=self.pbc)

        return ang

    def pb_br_pb_angle_hostrelative(self, ind, ax, coup_dir):
        orig_ang = self.pb_br_pb_angle(ind, ax)
        pln = (ax, coup_dir)
        pln_dist = self.br_distance_from_plane(ind, pln, ax)
        sign = np.sign(pln_dist)

        if sign == 1:
            ang = 360-orig_ang
        else:
            ang = orig_ang
        
        return ang

    def br_environment_coords(self, ind):
        o_host_frg = self.get_motif(ind).host
        
        ax_positions = np.array([
                o_host_frg.X[0],
                self.get_motif(self.shift_index(ind, 0, -1)).host.X[0],
                o_host_frg.X[1],
                self.get_motif(self.shift_index(ind, 1, -1)).host.X[1],
                o_host_frg.X[2],
                self.get_motif(self.shift_index(ind, 2, -1)).host.X[2],
            ])

        return ax_positions

    def br_distance_from_plane(self, ind, pln, br_dir):
        o_host_frg = self.get_motif(ind).host
        o_pos = o_host_frg.B
        pb_coords = self.pb_environment_coords(ind)
        br_coords = self.br_environment_coords(ind)

        all_coords = np.concatenate([pb_coords, br_coords])
        if self.pbc:
            all_coords = pbc_vectors(self.cell, o_pos, all_coords, relative=False)

        br_ind = 3+(2*br_dir)
        dist_dat = distance_of_point_on_plane(
                o_pos, all_coords[pln[0]], all_coords[pln[1]], all_coords[br_ind]
        )
        (ax_sign,) = {0,1,2} - set(pln)
        dist = dist_dat[0]*np.sign(dist_dat[1][ax_sign])

        return dist

def pbc_vectors(cell, v1, v2, orig_ind=0, relative=True):
    """Obtain shortest vectors with PBC for the given vectors"""

    pmglattice = PmgLattice(cell)
    v1_fr = pmglattice.get_fractional_coords(v1)
    v2_fr = pmglattice.get_fractional_coords(v2)
    
    # The first vector in v1 is taken as origin
    svecs = pbc_shortest_vectors(pmglattice, v1_fr, v2_fr)[orig_ind]
    if not relative:
        svecs = svecs+v1

    return svecs

def _dihedral_pbc(cell, p1, p2, p3, p4):
    """Compute a PBC-aware dihedral angle in degrees."""
    v1 = pbc_vectors(cell, p1, p2)[0]
    v2 = pbc_vectors(cell, p2, p3)[0]
    v3 = pbc_vectors(cell, p3, p4)[0]
    p23 = np.cross(v2, v3)
    p12 = np.cross(v1, v2)
    return np.degrees(
        np.arctan2(
            np.linalg.norm(v2) * np.dot(v1, p23),
            np.dot(p12, p23),
        )
    )

class LatticeTrajectory:
    """A trajectory of lattice objects for each structure in the original trajectory"""

    def __init__(self, lattice_list):
        self.lattices = lattice_list
        self.nt = len(lattice_list)
        
        init_lat = lattice_list[0]
        self.nx = init_lat.nx
        self.ny = init_lat.ny
        self.nz = init_lat.nz
        self.pbc = init_lat.pbc
    
    def compute_for_all(self, fn, *args, **kwargs):
        vals = []
        for lat in self.lattices:
            val = fn(lat, *args, **kwargs)
            vals.append(val)

        return vals
    
    def compute_for_all_cells(self, fn, *args, dtype=float, **kwargs):
        vals = []
        for lat in self.lattices:
            val = lat.compute_for_all(fn, *args, dtype=dtype, **kwargs)
            vals.append(val)

        return vals

    def compute_for_all_with_sys(self, sys_list, fn, *args, **kwargs):
        """Compute for all lattices with corresponding systems"""
        vals = []
        for lat, sys in zip(self.lattices, sys_list):
            val = fn(lat, sys, *args, **kwargs)
            vals.append(val)

        return vals

    def compute(self, func, **kwargs):
        """Apply a unit-cell function to every lattice in the trajectory."""
        from guesthost.compute import compute_trajectory

        return compute_trajectory(self, func, **kwargs)

    def compute_layerwise(self, func, **kwargs):
        """Apply a layerwise/global function to every lattice in the trajectory."""
        from guesthost.layerwise import compute_layerwise

        return compute_layerwise(self, func, **kwargs)

    def ma_self_correlation_P1(self):
        """MA self correlation with first order legendre polynomial P1"""

        fn_get_cn = lambda x, ind: x.get_motif(ind).ma.CN
        all_corr = []
        for ix, iy, iz in np.ndindex(self.nx, self.ny, self.nz):
            all_cn = [fn_get_cn(self.lattices[t], (ix, iy, iz)) for t in range(self.nt)]

            corr = self_correlation_P1(all_cn)
            all_corr.append(corr)

        return np.array(all_corr)

    def ma_self_correlation_P2(self):
        """MA self correlation with first order legendre polynomial P2"""

        fn_get_cn = lambda x, ind: x.get_motif(ind).ma.CN
        all_corr = []
        for ix, iy, iz in np.ndindex(self.nx, self.ny, self.nz):
            all_cn = [fn_get_cn(self.lattices[t], (ix, iy, iz)) for t in range(self.nt)]

            corr = self_correlation_P2(all_cn)
            all_corr.append(corr)

        return np.array(all_corr)
