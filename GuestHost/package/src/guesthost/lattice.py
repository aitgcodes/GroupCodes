import numpy as np
from guesthost.analysis.functions_modules import (
        proj2plane, self_correlation_P1, self_correlation_P2,
        distance_of_point_on_plane
)
from pymatgen.util.coord_cython import pbc_shortest_vectors
from pymatgen.core import Lattice as PmgLattice

class Motif:
    """Class representing one unit of a repeating collection of fragments"""
    
    def __init__(self, frg_list):
        self.fragments = frg_list
        self.nfrg = len(frg_list)

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


    def ma_torsion(self, ind, frg_ref, htyp="N"):
        ix, iy, iz = ind
        frg = self.motif_grid[ix, iy, iz].fragments[0]
        if htyp == "N":
            v1 = frg_ref.HN
        elif htyp == "C":
            v1 = frg_ref.HC

        val = frg.caltorsion(v1, htyp=htyp)

        return val

    def all_ma_torsions(self, lattice_ref, htyp="N"):
        torsions = np.zeros((self.nx, self.ny, self.nz))
        for ix, iy, iz in np.ndindex(self.nx, self.ny, self.nz):
            frg_ref = lattice_ref.motif_grid[ix,iy,iz].fragments[0]
            torsion = self.ma_torsion((ix, iy, iz), frg_ref, htyp=htyp)

            torsions[ix,iy,iz] = torsion

        return torsions

    def ma_orientation(self, ind, ax):
        "MA orientations for unitcell specified by `ind` with respect to axis, `ax`"""

        ix, iy, iz = ind
        cn_vec = self.motif_grid[ix, iy, iz].fragments[0].CN
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

    def get_motif(self, ind):
        ix, iy, iz = ind

        return self.motif_grid[ix, iy, iz]

    def pb_environment_coords(self, ind):
        ax_positions = np.array([
                self.get_motif(self.shift_index(ind, 0)).fragments[1].B,
                self.get_motif(self.shift_index(ind, 1)).fragments[1].B,
                self.get_motif(self.shift_index(ind, 2)).fragments[1].B
            ])

        return ax_positions

    def pb_environment_cell(self, ind):
        ix, iy, iz = ind
        o_host_frg = self.get_motif(ind).fragments[1]
        o_pos = o_host_frg.B
        
        ax_positions = self.pb_environment_coords(ind)

        if self.pbc:
            lcell = pbc_vectors(self.cell, o_pos, ax_positions).T
        else:
            lcell = (ax_positions - o_pos).T

        return lcell

    def pb_br_pb_angle(self, ind, ax):
        ix, iy, iz = ind
        hostfrg_orig = self.get_motif(ind).fragments[1]
        hostfrg_ax = self.get_motif(self.shift_index(ind, ax)).fragments[1]

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
        o_host_frg = self.get_motif(ind).fragments[1]
        
        ax_positions = np.array([
                o_host_frg.X[0],
                self.get_motif(self.shift_index(ind, 0, -1)).fragments[1].X[0],
                o_host_frg.X[1],
                self.get_motif(self.shift_index(ind, 1, -1)).fragments[1].X[1],
                o_host_frg.X[2],
                self.get_motif(self.shift_index(ind, 2, -1)).fragments[1].X[2],
            ])

        return ax_positions

    def br_distance_from_plane(self, ind, pln, br_dir):
        o_host_frg = self.get_motif(ind).fragments[1]
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

    def ma_self_correlation_P1(self):
        """MA self correlation with first order legendre polynomial P1"""

        fn_get_cn = lambda x, ind: x.motif_grid[ind[0],ind[1],ind[2]].fragments[0].CN
        all_corr = []
        for ix, iy, iz in np.ndindex(self.nx, self.ny, self.nz):
            all_cn = [fn_get_cn(self.lattices[t], (ix, iy, iz)) for t in range(self.nt)]

            corr = self_correlation_P1(all_cn)
            all_corr.append(corr)

        return np.array(all_corr)

    def ma_self_correlation_P2(self):
        """MA self correlation with first order legendre polynomial P2"""

        fn_get_cn = lambda x, ind: x.motif_grid[ind[0],ind[1],ind[2]].fragments[0].CN
        all_corr = []
        for ix, iy, iz in np.ndindex(self.nx, self.ny, self.nz):
            all_cn = [fn_get_cn(self.lattices[t], (ix, iy, iz)) for t in range(self.nt)]

            corr = self_correlation_P2(all_cn)
            all_corr.append(corr)

        return np.array(all_corr)
