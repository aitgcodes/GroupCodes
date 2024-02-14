from ase.io import read, write
import numpy as np
import sys
#from orderparams import listplotcorr as lpc

class MethylAmmonium:

    def __init__(self,atlst):

        self.nat = len(atlst)
        self.atlst = atlst
        self.C = None
        self.HC = None
        self.N = None
        self.HN = None
        self.CN = None

    def assign(self,pos):

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
        else:
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

        th = th/float(n)

        return th

class Host:

    def __init__(self,atlst):

        self.nat = len(atlst)
        self.atlst = atlst
        self.B = None
        self.X = None

    def assign(self,pos):

        ### BX3 ion assigned as 0:B, 1-3:X (in-cell, x<y<z)
        self.B = pos[self.atlst[0]]

        self.X = np.zeros((3,3),dtype='float64')
        for i in range(3):
            self.X[i] = pos[self.atlst[1+i]]

        self.Xo = np.zeros((3,3),dtype='float64') # out-cell X part of Octahedra (-x<-y<-z)

    def buildOh(self,hostlat):

        for ix in range(3):
            self.Xo[ix]=hostlat.bshft(ix).X[ix] 

        return 

    def OhDist(self):


#    def OhTilt:

class Trajectory:

    def __init__(self,fname):

        self.coords = []
        self.sym = []
        self.readxyz(fname)
        self.nat = len(self.coords[0])
        self.nt = len(self.coords)

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

    def fragmentize(self,ncells,atlst,fragment):

        ndiv = len(atlst)
        ucell = self.nat//ncells

        molt = []

        for atoms in self.coords:
            molecule = [fragment(atlst) for i in range(ncells)]
            for icell in range(ncells):
                molecule[icell].assign(atoms[icell*ucell:(icell+1)*ucell]) 
            molt.append(molecule)

        return molt
 
