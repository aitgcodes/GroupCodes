from gpaw.mpi import world
from gpaw import GPAW
from ase.io import write
from gpaw import restart
import numpy as np
from ase.units import Bohr
import os
import sys
import subprocess
sys.path.insert(0, os.path.dirname(__file__))

def ks_weights(wi_full_fpath, wa_full_fpath, wi_np_fpath, wa_np_fpath, i_p, a_p):

    imin = np.min(i_p)
    imax = np.max(i_p)
    amin = np.min(a_p)
    amax = np.max(a_p)
    print("imin=",imin)
    print("imax=",imax)
    print("amin=",amin)
    print("amax=",amax)

    data_wi_full_fpath = open(wi_full_fpath, "r")
    data_wa_full_fpath = open(wa_full_fpath, "r")
    data_wi_np_fpath = open(wi_np_fpath, "r")
    data_wa_np_fpath = open(wa_np_fpath, "r")

    wi_full = np.zeros(imax-imin+1)
    wa_full = np.zeros(amax-amin+1)
    wi_np = np.zeros(imax-imin+1)
    wa_np = np.zeros(amax-amin+1)

    wi_full_lines = data_wi_full_fpath.readlines()
    wa_full_lines = data_wa_full_fpath.readlines()
    wi_np_lines = data_wi_np_fpath.readlines()
    wa_np_lines = data_wa_np_fpath.readlines()

    # Read weights and normalize them
    for i in range(imin, imax+1):
        line_full = wi_full_lines[i-imin].strip().split()
        line_np = wi_np_lines[i-imin].strip().split()
        occw_full = float(line_full[2])
        occw_np = float(line_np[2])
        wi_full[i-imin] = occw_full / occw_full
        wi_np[i-imin] = occw_np / occw_full

    for a in range(amin, amax+1):
        line_full = wa_full_lines[a-amin].strip().split()
        line_np = wa_np_lines[a-amin].strip().split()
        unoccw_full = float(line_full[2])
        unoccw_np = float(line_np[2])
        wa_full[a-amin] = unoccw_full / unoccw_full 
        wa_np[a-amin] = unoccw_np / unoccw_full

    # Wave function weights for the occupied and unoccupied states in left (L) and right (R) regions of space
    wi_L = np.abs(wi_np)
    wi_R = np.abs(np.subtract(1, wi_L))
    wa_L = np.abs(wa_np)
    wa_R = np.abs(np.subtract(1, wa_L))

    # write the normalized weights of L and R regions in a file
    fout_i = open("wi_normalized.dat", "w")
    fout_a = open("wa_normalized.dat", "w")

    for i in range(imin, imax+1):
        fout_i.write("State"+"  "+str(i)+"  "+str(wi_L[i-imin]+wi_R[i-imin])+"  "+str(wi_L[i-imin])+"  "+str(wi_R[i-imin]))
        fout_i.write("\n")
    for a in range(amin, amax+1):
        fout_a.write("State"+"  "+str(a)+"  "+str(wa_L[a-amin]+wa_R[a-amin])+"  "+str(wa_L[a-amin])+"  "+str(wa_R[a-amin]))
        fout_a.write("\n")

    return [wi_L,wi_R,wa_L,wa_R]
