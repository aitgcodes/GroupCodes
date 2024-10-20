import os
from ase.io import write
from ase.units import Hartree, Bohr, AUT
from ase.io import Trajectory
from ase.parallel import parprint

name = 'opt'

traj = Trajectory(name + '.traj', 'r')

nsteps = len(traj)

print("nsteps =", nsteps)

string = 'structure'

#get current working directory and make a scratch 
#directory
path = os.getcwd()
path = path + '/scratch'
if not os.path.exists(path): os.makedirs(path)

#output file name
outFileName = 'trajectory.xyz'
#write each structure from the .traj file in .xyz format
for i in range(nsteps):
    atoms = traj[i]
    string = 'structure%03d' % (i,) +'.xyz'
    outStruct = os.path.join(path, string)
    write(outStruct, atoms)
#combines all structures in one trajectory 
#file
    inFile = open(os.path.join(path, 'structure%03d' %
                  (i,)  +'.xyz'), 'r')
    fileStr = inFile.read()
    outFile = open(outFileName, 'a')
    outFile.write(fileStr)
