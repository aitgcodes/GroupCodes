import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from ase.units import Hartree

file_np_tot = "Ag_tot.dat"
file_np_s = "Ag_s.dat"
file_np_p = "Ag_p.dat"
file_np_d = "Ag_d.dat"

dat_np_tot = np.loadtxt(file_np_tot)
dat_np_s = np.loadtxt(file_np_s)
dat_np_p = np.loadtxt(file_np_p)
dat_np_d = np.loadtxt(file_np_d)

energy = dat_np_s[:,0]

dos_np_tot = dat_np_tot[:,1]
dos_np_s = dat_np_s[:,1]
dos_np_p = dat_np_p[:,1]
dos_np_d = dat_np_d[:,1]

ft=16

# plot
plt.figure(figsize=(6.0, 5.4))
plt.plot(energy, dos_np_tot, label='Ag (total)', lw=2)
plt.plot(energy, dos_np_s, label='Ag (s)', lw=2)
plt.plot(energy, dos_np_p, label='Ag (p)', lw=2)
plt.plot(energy, dos_np_d, label='Ag (d)', lw=2)

plt.axis(ymin=0., ymax=180, xmin=-6.0, xmax=4.0, )
plt.vlines(x=0., ymin=0,ymax=180, colors='k', ls=':', lw=2)
plt.xlabel(r'$E - E_F \ \rm{(eV)}$', fontsize=ft)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('DOS', fontsize=ft)
plt.legend(loc='best', bbox_to_anchor=(0.45, 0.5, 0.5, 0.5), fontsize=14)
plt.savefig('Ag55_dos.png')
plt.show()
