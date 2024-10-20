import os
import sys
import numpy as np

data_1 = np.loadtxt('Ag147_4_ACF.dat')
data_2 = np.loadtxt('Ag55_4_ACF.dat')
data_3 = np.loadtxt('Au147_4_ACF.dat')
data_4 = np.loadtxt('Au55_4_ACF.dat')

Z_se = 34
Z_mo = 42
Z_ag = 47
Z_au = 79

fou = open('output.dat', 'w')

# Total charge on independent nanoparticles and TMD flake
Q_tot_ag147 = 147*Z_ag
Q_tot_ag55 = 55*Z_ag
Q_tot_au147 = 147*Z_au
Q_tot_au55 = 55*Z_au
Q_tot_tmd = 36*Z_se + 10*Z_mo

# Charge paritioning in Ag147_TMD
Q_ag147 = np.sum(data_1[0:147, 4])
Q_tmd = np.sum(data_1[147:193, 4])
dQ_ag147 = Q_tot_ag147 - Q_ag147
dQ_tmd   = Q_tot_tmd - Q_tmd
fou.write("# Charge paritioning in Ag147_TMD\n")
fou.write("dQ_ag147 = "+"   "+str(dQ_ag147))
fou.write("\n")
fou.write("dQ_tmd = "+"   "+str(dQ_tmd))
fou.write("\n")
# Charge paritioning in Ag55_TMD
Q_ag55 = np.sum(data_2[0:55, 4])
Q_tmd = np.sum(data_2[55:101, 4])
dQ_ag55 = Q_tot_ag55 - Q_ag55
dQ_tmd   = Q_tot_tmd - Q_tmd
fou.write("# Charge paritioning in Ag55_TMD\n")
fou.write("dQ_ag55 = "+"   "+str(dQ_ag55))
fou.write("\n")
fou.write("dQ_tmd = "+"   "+str(dQ_tmd))
fou.write("\n")
# Charge paritioning in Au147_TMD
Q_au147 = np.sum(data_3[0:147, 4])
Q_tmd = np.sum(data_3[147:193, 4])
dQ_au147 = Q_tot_au147 - Q_au147
dQ_tmd   = Q_tot_tmd - Q_tmd
fou.write("# Charge paritioning in Au147_TMD\n")
fou.write("dQ_au147 = "+"   "+str(dQ_au147))
fou.write("\n")
fou.write("dQ_tmd = "+"   "+str(dQ_tmd))
fou.write("\n")
# Charge paritioning in Au55_TMD
Q_au55 = np.sum(data_4[0:55, 4])
Q_tmd = np.sum(data_4[55:101, 4])
dQ_au55 = Q_tot_au55 - Q_au55
dQ_tmd   = Q_tot_tmd - Q_tmd
fou.write("# Charge paritioning in Au55_TMD\n")
fou.write("dQ_au55 = "+"   "+str(dQ_au55))
fou.write("\n")
fou.write("dQ_tmd = "+"   "+str(dQ_tmd))
fou.write("\n")
