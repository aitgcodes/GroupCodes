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

def ks_weights(occ_weight_fpath, unocc_weight_fpath, imin, imax, amin, amax):

    f1 = open(occ_weight_fpath, "r")
    f2 = open(unocc_weight_fpath, "r")

    wi = np.zeros(imax-imin+1)
    wa = np.zeros(amax-amin+1)

    occw_lines = f1.readlines()
    unoccw_lines = f2.readlines()

    print("imin=",imin)
    print("imax=",imax)
    print("amin=",amin)
    print("amax=",amax)

    for i in range(imin, imax+1):
        line = occw_lines[i-imin].strip().split()
        occw = float(line[2])
        wi[i-imin] = occw

    for a in range(amin, amax+1):
        line = unoccw_lines[a-amin].strip().split()
        unoccw = float(line[2])
        wa[a-amin] = unoccw

    # Wave function weights for the occupied and unoccupied states in left (L) and right (R) regions of space
    wi_L = np.abs(wi)
    wi_R = np.abs(np.subtract(1, wi_L))
    wa_L = np.abs(wa)
    wa_R = np.abs(np.subtract(1, wa_L))

    # normalize weights
    wi_L = wi_L/(wi_L+wi_R)
    wi_R = wi_R/(wi_L+wi_R)
    wa_L = wa_L/(wa_L+wa_R)
    wa_R = wa_R/(wa_L+wa_R)
    

    f1.close()
    f2.close()
    return [wi_L,wi_R,wa_L,wa_R]
