import os
import sys
from pathlib import Path
import numpy as np
from gpaw.tddft.units import au_to_fs


def extract_wtpdist(wtpdist_fpath, out_fpath):
    wtpdist = np.load(wtpdist_fpath)
    file_output = open(out_fpath, "w")

    time_t = wtpdist["time_t"]
    energy_o = wtpdist["energy_o"]
    energy_u = wtpdist["energy_u"]

    Ph_LR_to = wtpdist["Ph_LR_to"]
    Ph_LL_to = wtpdist["Ph_LL_to"]
    Ph_RR_to = wtpdist["Ph_RR_to"]
    Ph_RL_to = wtpdist["Ph_RL_to"]
    
    Pe_LR_tu = wtpdist["Pe_LR_tu"]
    Pe_LL_tu = wtpdist["Pe_LL_tu"]
    Pe_RR_tu = wtpdist["Pe_RR_tu"]
    Pe_RL_tu = wtpdist["Pe_RL_tu"]

    file_output.write("# Weighted transition probability distribution\n")
    file_output.write("# energy_o    energy_u  Ph_LR_to   Ph_LL_to  Ph_RR_to  Ph_RL_to   Pe_LR_tu   Pe_LL_tu   Pe_RR_tu   Pe_RL_tu ")
    file_output.write("\n")

    n_eo = len(energy_o)

    for i in range(n_eo):
        file_output.write(str(energy_o[i])+"     "+str(energy_u[i])+"     "+str(Ph_LR_to[0, i])+"     "+str(Ph_LL_to[0, i])+"     "+str(Ph_RR_to[0, i])
                +"     "+str(Ph_RL_to[0, i])+"     "+str(Pe_LR_tu[0, i])+"     "+str(Pe_LL_tu[0, i])+"     "+str(Pe_RR_tu[0, i])+"     "+str(Pe_RL_tu[0, i]))
        file_output.write("\n")
    file_output.close()

if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType, IntOrStrType

    parser = argparse.ArgumentParser()
    parser.add_argument('wtpdist_fpath', type=FilePathType)
    parser.add_argument('out_fpath', type=FilePathType)

    args = parser.parse_args()

    extract_wtpdist(args.wtpdist_fpath, args.out_fpath)
