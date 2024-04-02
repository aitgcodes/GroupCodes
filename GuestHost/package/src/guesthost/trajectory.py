from ase.io import read, write
import copy
import numpy as np
import sys
from guesthost.lattice import Motif, Lattice

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

    def assign(self,pos,indices=None):

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

    def planify(self,a):

        rmid = np.mean(a,axis=0)
        b = a - rmid

        plv = np.cross(b[0,:],b[1,:])
        plv = plv/np.linalg.norm(plv)

        return (b,plv)

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

class Host:

    def __init__(self,atlst):

        self.nat = len(atlst)
        self.atlst = atlst
        self.indices = None
        self.B = None
        self.X = None

    def assign(self,pos, indices=None):

        ### BX3 ion assigned as 0:B, 1-3:X (in-cell, x<y<z)
        self.B = pos[self.atlst[0]]

        self.X = np.zeros((3,3),dtype='float64')
        for i in range(3):
            self.X[i] = pos[self.atlst[1+i]]

        self.Xo = np.zeros((3,3),dtype='float64') # out-cell X part of Octahedra (-x<-y<-z)

        self.indices = indices[self.atlst]

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

    def readase(self, fname):
        # Read using ase
        atoms_list = read(fname, index=":")

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

        for atoms in self.coords:
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
                    molecule[icell].assign(atoms[idx], indices=np.array(idx))
            elif ordering == "unitcell":
                for icell in range(ncells):
                    idx = list(range(icell*ucell,(icell+1)*ucell))
                    molecule[icell].assign(atoms[idx], indices=np.array(idx))
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
