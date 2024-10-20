import os
import numpy as np
from ase.io import write
from gpaw import GPAW
from gpaw.tddft.units import au_to_eV
from gpaw.lcaotddft import LCAOTDDFT
from gpaw.lcaotddft.densitymatrix import DensityMatrix
from gpaw.lcaotddft.frequencydensitymatrix import FrequencyDensityMatrix
#from gpaw.lcaotddft.ksdecomposition import KohnShamDecomposition
from ksdecomposition import KohnShamDecomposition
import os

cubepath="cubes"
if not os.path.exists(cubepath):
    os.makedirs(cubepath)

def read_rho(time, tag):
    fpath = os.path.join("./rho", 't%09.1f%s.npy' % (time, tag))
    if not os.path.exists(fpath):
            raise RuntimeError('File missing: %s' % fpath)
    rho_p = np.load(fpath)
    return rho_p

calc= GPAW('unocc.gpw', txt=None)
calc.initialize_positions()  # Initialize in order to calculate density
dmat = DensityMatrix(calc)
ksd = KohnShamDecomposition(calc, 'ksd.ulm')

rho_MM_0 = read_rho(0, "")

list_t_0 = [i for i in range(11400, 11600, 20)]
list_t_1 = [i for i in range(17300, 17500, 20)]
list_t_2 = [25000]
#
list_t = [list_t_0, list_t_1, list_t_2]

for i in range(len(list_t)):
    for t in list_t[i]:
        rho_MM_t = read_rho(t, "")
        rho_up = rho_MM_t - rho_MM_0
        rho_g = ksd.get_density(calc.wfs, rho_up.imag)
        t_fs = t/1000
        write(f'{cubepath}/ind_{t_fs:.2f}fs.cube', calc.atoms, data=rho_g)
