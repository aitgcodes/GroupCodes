from ase.io import read, write
import itertools
import copy
import numpy as np
import sys
from guesthost.lattice import Motif, Lattice
from guesthost.unitcells import (
    get_unitcell_indexdata,
    host_indices_from_unitcell_data,
    ma_indices_from_unitcell_data,
    unit_order_from_reference,
    unit_order_from_unitcell_data,
    unit_orders_from_unitcell_data,
)

class MethylAmmonium:

    def __init__(self,atlst):
        """Class representing a MethylAmmonium molecule in a trajectory

        Attributes
        ---------
        nat : Number of atoms
        atlst : Indices of atoms in the molecule

        """

        self.nat = len(atlst)
        self.atlst = atlst
        self.indices = None
        self.C = None
        self.HC = None
        self.N = None
        self.HN = None
        self.CN = None

    def assign(self,pos,indices=None,atoms=None):

        ### MA assigned as 0:C, 1-3:H (of C), 4:N, 5-7:H (of N)
        self.C = pos[self.atlst[0]]

        self.HC = np.zeros((3,3),dtype='float64')
        for i in range(3):
            self.HC[i] = pos[self.atlst[1+i]]

        self.N = pos[self.atlst[4]]

        self.HN = np.zeros((3,3),dtype='float64')
        for i in range(3):
            self.HN[i] = pos[self.atlst[5+i]]

        d = self.N - self.C
        self.CN = d/np.linalg.norm(d)
        
        self.indices = indices[self.atlst]
        self.atoms = atoms[self.atlst]

    def planify(self,a):

        rmid = np.mean(a,axis=0)
        b = a - rmid

        plv = np.cross(b[0,:],b[1,:])
        plv = plv/np.linalg.norm(plv)

        return (b,plv)

    def costheta(self,a,b):
   
        a1 = a/np.linalg.norm(a)
        b1 = b/np.linalg.norm(b)
        cst = np.dot(a1,b1)
        #print(cst)
        if cst < 0:
           cst = max(cst,-1.0)
        elif cst > 0:
           cst = min(cst,1.0)

        return cst

    def angle(self,a,b):
   
        a1 = a/np.linalg.norm(a)
        b1 = b/np.linalg.norm(b)
        cst = np.dot(a1,b1)
        #print(cst)
        if cst < 0:
           cst = max(cst,-1.0)
        elif cst > 0:
           cst = min(cst,1.0)

        #print(cst)
        theta = np.degrees(np.arccos(cst))

        return theta
    
    def caltorsion(self,v1,htyp="N"):
### Expects two arrays each of shape (3x3) with the coordinates of the 3 hydrogens

        if htyp == "N":
           v2 = self.HN
        elif htyp == "C":
           v2 = self.HC

        n = len(v2)

        (v1int,pln1) = self.planify(v1)
        (v2int,pln2) = self.planify(v2)

        diff = v2int-v1int

    ### Cos of angle between vectors at two different times
        th = 0.0

        for i in range(n):
            vec = np.cross(diff[i,:],v1int[i,:])
            norm = max(np.linalg.norm(vec),1.e-8)
            #print(norm)
            
            ph = np.dot(vec,pln1)/norm
            ph = ph/max(abs(ph),1.e-8)
            th += np.round(ph)*self.angle(v1int[i,:],v2int[i,:])
            # th += self.angle(v1int[i,:],v2int[i,:])

        th = th/float(n)

        return th

    def caltorsion_atan(self, v1, htyp="N"):
        """
        Compute the signed average torsional rotation between two snapshots of an MA end-group.
        Uses θ_i = atan2( n̂ · (u_i × v_i), u_i · v_i ) for each H and averages over the triad.
        
        Args:
            v1  : (3,3) array of H coordinates at time t1 (rows = H1,H2,H3)
            htyp: "N" -> compare to self.HN; "C" -> compare to self.HC (coordinates at time t2)

        Returns:
            Mean signed angle (radians) in (-pi, pi].
        """
        # Choose the comparison triad (time t2)
        if htyp == "N":
            v2 = self.HN
        elif htyp == "C":
            v2 = self.HC
        else:
            raise ValueError("htyp must be 'N' or 'C'")

        # Center each triad and get the initial plane normal n̂ from v1
        v1int, n_hat = self.planify(v1)
        v2int, _     = self.planify(v2)  # plane of v2 not needed for the sign

        n = v1int.shape[0]
        eps = 1e-12
        theta_sum = 0.0

        for i in range(n):
            u = v1int[i, :]
            v = v2int[i, :]

            # Normalize to remove any scale dependence
            u /= max(np.linalg.norm(u), eps)
            v /= max(np.linalg.norm(v), eps)

            sin_term = np.dot(n_hat, np.cross(u, v))  # oriented sine component
            cos_term = np.dot(u, v)                   # cosine component

            theta_i = np.degrees(np.arctan2(sin_term, cos_term))  # signed angle in (-180, 180]
            theta_sum += theta_i

        return theta_sum / float(n)

    def caltorsion_removetilts(self, v1, htyp="N", match_rows=False, return_tilt=False):
        """
        Signed torsion of a 3-H triad about its local C3 axis, cleaned of tumbling.
        Implements: project the later vector triad into the initial plane, then
        compute the oriented in-plane rotation and average over the three H's.

        Args
        ----
        v1 : (3,3) array
            H coordinates at time t1 (rows = H1,H2,H3).
        htyp : {'N','C'}
            Compare to self.HN (default) or self.HC at time t2.
        match_rows : bool
            If True, permute v2 rows to best match v1 (helps if H ordering changes).
        return_tilt : bool
            If True, also return the tilt angle between the two planes.

        Returns
        -------
        theta : float
            Mean signed torsion in radians in (-pi, pi].
        (tilt) : float, optional
            Tilt between the plane normals (0..pi).
        """
        # --- choose comparison triad (time t2)
        if   htyp == "N": v2 = self.HN
        elif htyp == "C": v2 = self.HC
        else: raise ValueError("htyp must be 'N' or 'C'")

        # --- center and get plane normals
        v1int, n1 = self.planify(v1)   # initial triad & its normal (sign sets convention)
        v2int, n2 = self.planify(v2)   # later triad (normal only used for 'tilt' output)

        eps = 1e-12

        # --- optionally fix hydrogen order by trying all 6 permutations
        if match_rows:
            best_score, best_perm = -np.inf, (0,1,2)
            u_norms = np.linalg.norm(v1int, axis=1) + eps
            for perm in itertools.permutations((0,1,2)):
                w = v2int[list(perm)]
                w_norms = np.linalg.norm(w, axis=1) + eps
                # similarity score = sum of |cosine| between paired vectors
                score = np.sum(np.abs(np.sum(v1int*w, axis=1) / (u_norms*w_norms)))
                if score > best_score:
                    best_score, best_perm = score, perm
            v2int = v2int[list(best_perm)]

        # --- project later vectors into the initial plane to remove tumbling
        # row-wise: v2_proj_i = v2_i - (v2_i·n1) n1
        proj_coeff = v2int @ n1              # shape (3,)
        v2proj = v2int - proj_coeff[:, None] * n1

        # --- compute oriented in-plane angle per H and average
        theta_sum = 0.0
        for u, v in zip(v1int, v2proj):
            u = u / (np.linalg.norm(u) + eps)
            v = v / (np.linalg.norm(v) + eps)
            sin_term = np.dot(n1, np.cross(u, v))   # oriented sine in the n1 direction
            cos_term = np.dot(u, v)                 # cosine
            theta_sum += np.degrees(np.arctan2(sin_term, cos_term))

        theta = theta_sum / 3.0

        if return_tilt:
            # tilt between plane normals (for diagnostics)
            sin_tilt = np.linalg.norm(np.cross(n1, n2))
            cos_tilt = np.dot(n1, n2)
            tilt = np.degrees(np.arctan2(sin_tilt, cos_tilt))
            return theta, tilt

        return theta

class Host:

    def __init__(self,atlst):

        self.nat = len(atlst)
        self.atlst = atlst
        self.indices = None
        self.B = None
        self.X = None

    def assign(self,pos, indices=None, atoms=None):

        ### BX3 ion assigned as 0:B, 1-3:X (in-cell, x<y<z)
        self.B = pos[self.atlst[0]]

        self.X = np.zeros((3,3),dtype='float64')
        for i in range(3):
            self.X[i] = pos[self.atlst[1+i]]

        self.Xo = np.zeros((3,3),dtype='float64') # out-cell X part of Octahedra (-x<-y<-z)

        self.indices = indices[self.atlst]
        self.atoms = atoms[self.atlst]

    def buildOh(self,hostlat):

        for ix in range(3):
            self.Xo[ix]=hostlat.bshft(ix).X[ix] 

        return 

#    def OhDist(self):


#    def OhTilt:

# Create indices to order supercell based on unitcell
# in case the supercell is ordered based on element types
def get_supercell_indices(unit_order, ncells=64):
    ordered_inds = []
    for i in range(ncells):
        c, h_c, n, h_n, pb, br = copy.deepcopy(unit_order)
        c += i
        h_c += 6*i
        n += i
        h_n += 6*i
        pb += i
        br += 3*i
        current_inds = [
                c, h_c[0], h_c[1], h_c[2],
                n, h_n[0], h_n[1], h_n[2],
                pb, br[0], br[1], br[2]
            ]
        ordered_inds += current_inds

    return ordered_inds

# Default order in unitcell
unit_order = [
        320,                        # C
        np.array([384, 385, 386]),  # H_C
        256,                        # N
        np.array([387, 388, 389]),  # H_N
        0,                          # Pb
        np.array([64, 65, 66])      # Br
    ]

default_ordered_inds = get_supercell_indices(unit_order)

class Trajectory:

    def __init__(self,fname, order=False):
        """A class representing a trajectory containing coordinates,
        cells, symbols, number of atoms and number of steps."""

        self.coords = []
        self.sym = []
        self.cells = []
        
        # self.readxyz(fname)
        self.readase(fname)

        self.nat = len(self.coords[0])
        self.nt = len(self.coords)

    @classmethod
    def from_atoms_list(cls, atoms_list):
        """Create a Trajectory from an in-memory sequence of ASE Atoms."""
        obj = cls.__new__(cls)
        obj.coords = []
        obj.sym = []
        obj.cells = []
        obj.atoms_list = list(atoms_list)
        for atoms in obj.atoms_list:
            obj.coords.append(atoms.positions)
            obj.cells.append(atoms.get_cell().array)
        obj.sym = np.array(obj.atoms_list[0].get_chemical_symbols())
        obj.nat = len(obj.coords[0])
        obj.nt = len(obj.coords)
        return obj

    def readase(self, fname):
        # Read using ase
        atoms_list = read(fname, index=":")
        self.atoms_list = atoms_list

        for atoms in atoms_list:
            self.coords.append(atoms.positions)
            self.cells.append(atoms.get_cell().array)

        self.sym = np.array(atoms_list[0].get_chemical_symbols())

    def readxyz(self,fname):
 
        fp=open(fname,'r')
        lines=fp.readlines()
        fp.close()
        nat = int(lines[0].strip())
        nlines = len(lines)
        nt = nlines//(nat+2)
        print(nt)

        for i in range(nt):
            ibeg = i*(nat+2)
            pos = np.zeros((nat,3),dtype='float64')
            sym = []
            for iat in range(nat):
                sym.append(lines[ibeg+iat+2].strip().split()[0])
                pos[iat,:] = np.array(lines[ibeg+iat+2].strip().split()[1:]).astype(float)

            self.coords.append(pos)

            if i == 0:
               self.sym = sym

            #print("Read line {}.".format(it))

    def subtraj(self,traj,atlst):

        strj = []
        nat = len(atlst)
        atpos = np.zeros((nat,3),dtype='float64')
        for molecule in traj.coords:
            for iat in len(nat):
                atpos[iat,:] = molecule[atlst[iat],:] 
    
            strj.append(atpos)

        return strj

    def fragmentize(self,ncells,atlst,fragment,ordering="unitcell",unit_order=None):
        """Returns list of fragment objects for each step in the trajectory"""

        ndiv = len(atlst)
        ucell = self.nat//ncells
        atlst = np.array(atlst)

        molt = []

        for i, atoms in enumerate(self.coords):
            molecule = [fragment(atlst) for i in range(ncells)]
            if ordering == "type":
                if unit_order is None:
                    raise TypeError("Expected `unit_order` specifying the indices of atoms in home unit cell")
                unit_order = np.array(unit_order)
                unit_syms = self.sym[unit_order]
                unqspec, spec_counts = np.unique(unit_syms, return_counts=True)
                count_dict = dict(zip(unqspec, spec_counts))
                for icell in range(ncells):
                    idx = []
                    for at_ind, at_sym in zip(unit_order, unit_syms):
                        at_ind_c = at_ind + (icell * count_dict[at_sym])
                        idx.append(at_ind_c)
                    molecule[icell].assign(atoms[idx], indices=np.array(idx), atoms=self.atoms_list[i][idx])
            elif ordering == "unitcell":
                for icell in range(ncells):
                    idx = list(range(icell*ucell,(icell+1)*ucell))
                    molecule[icell].assign(atoms[idx], indices=np.array(idx), atoms=self.atoms_list[i][idx])
            molt.append(molecule)

        return molt

    def create_lattice(
            self, ncells, atlst_frg, frgtype_list, supercell_size, lattice_type,
            ordering="unitcell", unit_order=None
        ):
        """Create Lattice object for the given fragment types

        Parameters
        ----------
        ncells: Number of unit cells
        atlst_frg: List of lists containing the indices of each fragment
            in the parent unit cell
        frgtype_list: Types of each fragment
        """
        nfrg = len(atlst_frg)
        nx, ny, nz = supercell_size

        allfrg_lists = []
        for atlst, fragment in zip(atlst_frg, frgtype_list):
            frg_lists = self.fragmentize(
                ncells, atlst, fragment, ordering=ordering, unit_order=unit_order
            )
            allfrg_lists.append(frg_lists)

        lattices = []
        for it in range(self.nt):
            allmotifs = []
            for icell in range(ncells):
                motif_frgs = []
                for ifrg in range(nfrg):
                    motif_frgs.append(allfrg_lists[ifrg][it][icell])

                motif = Motif(motif_frgs)
                allmotifs.append(motif)

            lattice = lattice_type(allmotifs, self.cells[it], nx, ny, nz)
            lattices.append(lattice)

        return lattices

    def create_hplattice_from_unitcell_data(self, unitcell_data, supercell_size=None):
        """Create HPLattice objects from explicit unit-cell index data."""
        from guesthost.lattice import HPLattice

        if supercell_size is None:
            ncells = len(unitcell_data)
            n = round(ncells ** (1 / 3))
            if n ** 3 != ncells:
                raise ValueError("Pass supercell_size=(nx, ny, nz) for non-cubic systems.")
            supercell_size = (n, n, n)
        else:
            ncells = int(np.prod(supercell_size))

        if len(unitcell_data) != ncells:
            raise ValueError(
                f"unitcell_data has {len(unitcell_data)} cells, expected {ncells}."
            )

        nx, ny, nz = supercell_size
        lattices = []
        local_ma_inds = np.arange(8)
        local_host_inds = np.arange(4)

        for iframe, frame in enumerate(self.atoms_list):
            motifs = []
            positions = self.coords[iframe]
            for udata in unitcell_data:
                ma_global = np.array(ma_indices_from_unitcell_data(udata), dtype=int)
                host_global = np.array(host_indices_from_unitcell_data(udata), dtype=int)

                ma = MethylAmmonium(local_ma_inds)
                ma.assign(positions[ma_global], indices=ma_global, atoms=frame[ma_global])

                host = Host(local_host_inds)
                host.assign(positions[host_global], indices=host_global, atoms=frame[host_global])

                br_cage_global = np.array(udata["br_inds"], dtype=int)
                motifs.append(
                    Motif(
                        [ma, host],
                        unitcell_data=udata,
                        br_cage_coords=positions[br_cage_global],
                        br_cage_indices=br_cage_global,
                    )
                )

            lattices.append(HPLattice(motifs, self.cells[iframe], nx, ny, nz))

        return lattices

    def create_hplattice_from_reference(self, reference=None, supercell_size=None, **kwargs):
        """Infer unit-cell index data from a reference frame and create HPLattices."""
        if reference is None:
            reference = self.atoms_list[0]
        unitcell_data = get_unitcell_indexdata(reference, **kwargs)
        return self.create_hplattice_from_unitcell_data(unitcell_data, supercell_size=supercell_size)

    def create_hplattice(self, supercell_size=None, unitcell_data=None, reference=None, **kwargs):
        """Create HPLattice objects using unit-cell data or a reference structure."""
        if unitcell_data is not None:
            return self.create_hplattice_from_unitcell_data(unitcell_data, supercell_size=supercell_size)
        return self.create_hplattice_from_reference(
            reference=reference,
            supercell_size=supercell_size,
            **kwargs,
        )
