import os
import sys
from pathlib import Path
import numpy as np
from gpaw.tddft.units import au_to_fs


def extract_energy(energy_fpath, out_fpath):
    energy = np.load(energy_fpath)
    file_output = open(out_fpath, "w")

    time_t = energy["time_t"]/1000   # To convert time from as to fs
    E_t = energy["E_t"]
    Ec_t = energy["Ec_t"]
    Eq_t = energy["Eq_t"]
    Ep_t = energy["Ep_t"]

    file_output.write("# Energy contributions from pulse convolution\n")
    file_output.write("# Time     Total                         Coulomb                    - terms with q                   - terms with p")
    file_output.write("\n")
    
    n_t = len(time_t)
    
    for i in range(n_t):
        file_output.write(str(time_t[i])+"     "+str(E_t[i])+"     "+str(Ec_t[i])+"     "+str(Eq_t[i])+"     "+str(Ep_t[i]))
        file_output.write("\n")
    
    file_output.close()

if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType, IntOrStrType

    parser = argparse.ArgumentParser()
    parser.add_argument('energy_fpath', type=FilePathType)
    parser.add_argument('out_fpath', type=FilePathType)

    args = parser.parse_args()

    extract_energy(args.energy_fpath, args.out_fpath)
