# web-page: spec.png
import numpy as np
import matplotlib.pyplot as plt

data_dm = np.loadtxt('dmpulse.dat')
data_pulse = np.loadtxt('pulse.dat')

plt.figure(figsize=(6.0, 4.0), dpi=300)

plt.plot(0.0242*data_dm[:, 0], (data_dm[:, 4] - data_dm[0, 4]), label='Dipole moment oscillations')
plt.plot(0.0242*data_pulse[:, 0], 0.06*1.0e5*data_pulse[:, 1], label='Laser Pulse')
plt.xlabel('Time (fs)')
#plt.ylabel('Dipole moment (au)')
plt.xlim(0, 60)
#plt.ylim(ymin=1.0, ymax=2.4)
plt.legend(loc='best', bbox_to_anchor=(0.4, 0.5, 0.5, 0.5))
plt.tick_params(left = False, labelleft = False) 
plt.savefig('ag55_pfoa_dm_pulse.png')
plt.show()
