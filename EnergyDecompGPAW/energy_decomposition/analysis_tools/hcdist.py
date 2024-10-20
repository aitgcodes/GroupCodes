import os
import sys
from pathlib import Path
import numpy as np
from gpaw.tddft.units import au_to_fs


def extract_hcdist(hcdist_fpath, out_fpath):
    hcdist = np.load(hcdist_fpath)
    file_output = open(out_fpath, "w")

    time_t = hcdist["time_t"]
    energy_o = hcdist["energy_o"]
    energy_u = hcdist["energy_u"]

    dist_to = hcdist["dist_to"]
    dist_tu = hcdist["dist_tu"]

    file_output.write("# Hot carrier distribution\n")
    file_output.write("# energy_o      hole distribution      energy_u      electron distribution")
    file_output.write("\n")

    n_eo = len(energy_o)

    for i in range(n_eo):
        file_output.write(str(energy_o[i])+"     "+str(dist_to[0, i])+"     "+str(energy_u[i])+"     "+str(dist_tu[0, i]))
        file_output.write("\n")
    file_output.close()

if __name__ == '__main__':
    import argparse
    from argparse_util import ExistingPathType, FilePathType, IntOrStrType

    parser = argparse.ArgumentParser()
    parser.add_argument('energy_fpath', type=FilePathType)
    parser.add_argument('out_fpath', type=FilePathType)

    args = parser.parse_args()

    extract_hcdist(args.energy_fpath, args.out_fpath)
