import sys
import numpy as np
from ase.io import read
from trajproc import Trajectory, MethylAmmonium
from orderparams import listplotcorr as lpc
from orderparams import gridcorrvec as gcv

#### Read in the parameters
n = 4 ### No. of unit cells in each direction
delt = 0.00120945 ### Time step in trajectory
fname = sys.argv[1] ### Trajectory file name (XYZ format)
atln = sys.argv[2].strip().split() ### List of indices for fragment in each unit cell

atlst = []
for atl in atln:
    atlst.append(int(atl))

### Read in the trajectory
traj = Trajectory(fname)
#traj = read(fname)
print("Read in trajectory of length {} with {} atoms at each step!".format(traj.nt,traj.nat))

nunits = n*n*n
nt = traj.nt
print(nt,nunits)

### Extract coordinates of fragment objects
### This is a list of lists of objects of type MethylAmmonium
mas_t = traj.fragmentize(nunits,atlst,MethylAmmonium)
mas0 = mas_t[0]
print("Extracted fragments!")

tvals = np.linspace(0,(nt-1)*delt,nt)
print(tvals.size)

tort = []
madipt = []
for mas in mas_t:
    torsion = []
    madip = []
    for iu in range(nunits):
        torsion.append(mas[iu].caltorsion(mas0[iu].HN,htyp="N"))
        madip.append(mas[iu].CN)

    torsion = np.array(torsion)
    torsion = np.reshape(np.array(torsion),(n,n,n))
    madip = np.array(madip)
    tort.append(torsion)
    madipt.append(madip)

(data,err) = lpc(tort,tvals,mode="corr",avg="all")
fp=open('data.dat','w')
for i in range(len(data)):
    fp.write("{} {}\n".format(data[i,0],data[i,1]))
fp.close()

### Calculate spatial correlation of CN vecs
dirn = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(0,1,1),(1,0,1),(2,0,0)]
shp = (n,n,n)
corrt = []
for idir in dirn:
    corrt.append(gcv(madipt,idir,shp))

print(corrt)
