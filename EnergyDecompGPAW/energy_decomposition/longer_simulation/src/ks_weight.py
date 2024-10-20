#This codes integrates cube files upto a bounadary specified along the x , y and z axis respectively
# And it calculates the weight of specified number of Kohn-Sham states

def ks_weights(cube_dpath, out_fpath, bands_min, bands_max, lx, ly, lz):
    from gpaw.mpi import world
    from ase.io import write
    import numpy as np
    import os
    import sys
    sys.path.insert(0, os.path.dirname(__file__))

    import subprocess

    f = open(out_fpath, "w")

    for i in range(bands_min, bands_max):
        out = subprocess.Popen(["python3","src/cubeintegrator.py",cube_dpath+"gs_"+str(i)+".cube", lx,ly,lz],stdout=subprocess.PIPE,universal_newlines=True)
        out.wait()  # wait until the script has finished
        stdout_data, stderr_data = out.communicate()
        f.writelines(["State  "+str(i)+"  ",stdout_data])
    f.close()



if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('cube_dpath', type=ExistingPathType)
    parser.add_argument('out_fpath', type=FilePathType)
    parser.add_argument('--bands_min', type=int, default=0)
    parser.add_argument('--bands_max', type=int, default=0)
    parser.add_argument('--lx', type=str, default='-')
    parser.add_argument('--ly', type=str, default='-')
    parser.add_argument('--lz', type=str, default='-')
    args = parser.parse_args()

    ks_weights(args.cube_dpath, args.out_fpath, args.bands_min, args.bands_max, args.lx, args.ly, args.lz)
