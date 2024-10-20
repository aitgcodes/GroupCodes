import numpy as np

'''
  This code extracts the total projected dos on NP cluster and PFAS molecule
  and also on atomic orbitals
'''

# Define the list of atomic orbitals for Ag atom
angular=['s', 'p', 'd']
infile_test = "data/Ag_0_s.dat"
dat_test = open(infile_test, "r")
lines_test = dat_test.readlines()
e = np.zeros(len(lines_test))
#--------- Ag------------------
dos_np_tot = np.zeros(len(lines_test))
out_f_np_tot = open("Ag55_tot.dat", "w")
for c in angular:
    dos_t = np.zeros(len(e))
    outfile_np = "Ag_{}.dat".format(c)
    for a in range(55):
        infile = "data/Ag_{}_{}.dat".format(a, c)
        dat = open(infile, "r")
        lines = dat.readlines()
        for i in range(len(lines)):
            line = lines[i].strip().split()
            e[i] = float(line[0])
            dos = float(line[1])
            dos_t[i] += dos
    dos_np_tot += dos_t  # Total dos of NP
    out_f = open(outfile_np, "w")
    for i in range(len(e)):
        out_f.write(str(e[i])+"    "+str(dos_t[i])+"\n")
for i in range(len(e)):
    out_f_np_tot.write(str(e[i])+"    "+str(dos_np_tot[i])+"\n")

#--------- TMD------------------
#------total dos-
dos_tmd_tot = np.zeros(len(lines_test))
out_f_tmd_tot = open("MoSe2_4_tot.dat", "w")
for c in angular:
    dos_t = np.zeros(len(e))
    for a in range(55, 101, 1):
        infile = "data/tmd_{}_{}.dat".format(a, c)
        dat = open(infile, "r")
        lines = dat.readlines()
        for i in range(len(lines)):
            line = lines[i].strip().split()
            e[i] = float(line[0])
            dos = float(line[1])
            dos_t[i] += dos
    dos_tmd_tot += dos_t  # Total dos of PFOA
for i in range(len(e)):
    out_f_tmd_tot.write(str(e[i])+"    "+str(dos_tmd_tot[i])+"\n")

#---atom and orbital projected dos of tmd---------
#-----Se---
for c in angular:
    dos_t = np.zeros(len(e))
    outfile_Se = "Se_{}.dat".format(c)
    for a in range(55, 91, 1):
        infile = "data/tmd_{}_{}.dat".format(a, c)
        dat = open(infile, "r")
        lines = dat.readlines()
        for i in range(len(lines)):
            line = lines[i].strip().split()
            e[i] = float(line[0])
            dos = float(line[1])
            dos_t[i] += dos
    out_Se = open(outfile_Se, "w")
    for i in range(len(e)):
        out_Se.write(str(e[i])+"    "+str(dos_t[i])+"\n")
#-----Mo---
for c in angular:
    dos_t = np.zeros(len(e))
    outfile_Mo = "Mo_{}.dat".format(c)
    for a in range(91, 101, 1):
        infile = "data/tmd_{}_{}.dat".format(a, c)
        dat = open(infile, "r")
        lines = dat.readlines()
        for i in range(len(lines)):
            line = lines[i].strip().split()
            e[i] = float(line[0])
            dos = float(line[1])
            dos_t[i] += dos
    out_Mo = open(outfile_Mo, "w")
    for i in range(len(e)):
        out_Mo.write(str(e[i])+"    "+str(dos_t[i])+"\n")
