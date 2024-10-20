from gpaw.lcaotddft import LCAOTDDFT
from gpaw.lcaotddft.densitymatrix import DensityMatrix
from gpaw.lcaotddft.frequencydensitymatrix import FrequencyDensityMatrix
from gpaw.tddft.folding import frequencies

# Read the ground-state file
td_calc = LCAOTDDFT('gs.gpw')

# Attach analysis tools
dmat = DensityMatrix(td_calc)
freqs = frequencies([0.41, 1.65, 2.0, 2.50, 2.91, 3.35, 3.77, 4.05, 4.35, 4.52, 4.93, 5.14, 5.46], 'Gauss', 0.07)
fdm = FrequencyDensityMatrix(td_calc, dmat, frequencies=freqs)

# Replay the propagation
td_calc.replay(name='wf.ulm', update='none')

# Store the density matrix
fdm.write('fdm.ulm')
