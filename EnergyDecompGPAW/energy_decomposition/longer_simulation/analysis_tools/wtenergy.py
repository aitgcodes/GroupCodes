import os
import sys
from pathlib import Path
import numpy as np
from gpaw.tddft.units import au_to_fs


def extract_wtenergy(wtenergy_fpath, out_fpath):
    wtenergy = np.load(wtenergy_fpath)
    file_output = open(out_fpath, "w")

    time_t = wtenergy["time_t"]/1000   # To convert time from as to fs
    E_LL_t = wtenergy["E_LL_t"]
    E_RR_t = wtenergy["E_RR_t"]
    E_LR_t = wtenergy["E_LR_t"]
    E_RL_t = wtenergy["E_RL_t"]

    file_output.write("# Weighted Energy contributions from pulse convolution\n")
    file_output.write("# Time  E_LL_t   E_RR_t  E_LR_t  E_RL_t ")
    file_output.write("\n")
    
    n_t = len(time_t)
    
    for i in range(n_t):
        file_output.write(str(time_t[i])+"  "+str(E_LL_t[i])+"  "+str(E_RR_t[i])+"  "+str(E_LR_t[i])+"  "+str(E_RL_t[i]))
        file_output.write("\n")
    
    file_output.close()

if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType, IntOrStrType

    parser = argparse.ArgumentParser()
    parser.add_argument('wtenergy_fpath', type=FilePathType)
    parser.add_argument('out_fpath', type=FilePathType)

    args = parser.parse_args()

    extract_wtenergy(args.wtenergy_fpath, args.out_fpath)
