import numpy as np
from guesthost.analysis.functions_modules import (
        proj2plane, self_correlation_P1, self_correlation_P2
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
        self.pbc = pbc

    def compute_for_all(self, fn, dtype=float, *args, **kwargs):
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

    def pb_environment_cell(self, ind):
        ix, iy, iz = ind
        nx = self.nx
        ny = self.ny
        nz = self.nz
        o_host_frg = self.motif_grid[ix, iy, iz].fragments[1]
        o_pos = o_host_frg.B
        
        host_frgs = [
            self.motif_grid[(ix+1)%nx,iy,iz].fragments[1],
            self.motif_grid[ix,(iy+1)%ny,iz].fragments[1],
            self.motif_grid[ix,iy,(iz+1)%nz].fragments[1]
        ]
        # ax_positions = np.array([
        #         self.motif_grid[(ix+1)%nx,iy,iz].fragments[1].B,
        #         self.motif_grid[ix,(iy+1)%ny,iz].fragments[1].B,
        #         self.motif_grid[ix,iy,(iz+1)%nz].fragments[1].B
        #     ])
        ax_positions = np.array([
            host_frgs[0].B,
            host_frgs[1].B,
            host_frgs[2].B
        ])
        if self.pbc:
            lcell = pbc_vectors(self.cell, o_pos, ax_positions).T
        else:
            lcell = (ax_positions - o_pos).T

        return lcell

def pbc_vectors(cell, v1, v2):
    """Obtain shortest vectors with PBC for the given vectors"""

    pmglattice = PmgLattice(cell)
    v1_fr = pmglattice.get_fractional_coords(v1)
    v2_fr = pmglattice.get_fractional_coords(v2)
    
    # The first vector in v1 is taken as origin
    svecs = pbc_shortest_vectors(pmglattice, v1_fr, v2_fr)[0]

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
