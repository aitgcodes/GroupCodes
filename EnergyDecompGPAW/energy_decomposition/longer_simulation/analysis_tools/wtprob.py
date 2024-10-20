import os
import sys
from pathlib import Path
import numpy as np
from gpaw.tddft.units import au_to_fs


def extract_wtprob(wtprob_fpath, out_fpath_1, out_fpath_2):
    wtp = np.load(wtprob_fpath)
    file_output_1 = open(out_fpath_1, "w")
    file_output_2 = open(out_fpath_2, "w")

    time_t = wtp["time_t"]/1000   # To convert time from as to fs
    N_t = wtp["N_t"]
    N_LR_t = wtp["N_LR_t"]
    N_LL_t = wtp["N_LL_t"]
    N_RR_t = wtp["N_RR_t"]
    N_RL_t = wtp["N_RL_t"]
    
    P_LR_t = wtp["P_LR_t"]
    P_LL_t = wtp["P_LL_t"]
    P_RR_t = wtp["P_RR_t"]
    P_RL_t = wtp["P_RL_t"]

    file_output_1.write("# Number of carriers generated\n")
    file_output_1.write("# Time            N_t                    N_LR_t                  N_LL_t                N_RR_t              N_RL_t")
    file_output_1.write("\n")
    
    n_t = len(time_t)
    
    for i in range(n_t):
        file_output_1.write(str(time_t[i])+"     "+str(N_t[i])+"     "+str(N_LR_t[i])+"     "+str(N_LL_t[i])+"     "+str(N_RR_t[i])+"      "+str(N_RL_t[i]))
        file_output_1.write("\n")

    file_output_2.write("# Total probability of a partial process\n")
    file_output_2.write("# Time             P_LR_t                  P_LL_t                P_RR_t              P_RL_t")
    file_output_2.write("\n")
        
    for i in range(n_t):
        file_output_2.write(str(time_t[i])+"     "+str(P_LR_t[i])+"     "+str(P_LL_t[i])+"     "+str(P_RR_t[i])+"      "+str(P_RL_t[i]))
        file_output_2.write("\n")
    
    file_output_1.close()
    file_output_2.close()

if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType, IntOrStrType

    parser = argparse.ArgumentParser()
    parser.add_argument('wtprob_fpath', type=FilePathType)
    parser.add_argument('out_fpath_1', type=FilePathType)
    parser.add_argument('out_fpath_2', type=FilePathType)

    args = parser.parse_args()

    extract_wtprob(args.wtprob_fpath, args.out_fpath_1, args.out_fpath_2)
