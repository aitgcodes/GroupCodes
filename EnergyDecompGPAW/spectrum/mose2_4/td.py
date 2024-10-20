from ase.parallel import paropen
from gpaw.mpi import world
from gpaw.lcaotddft import LCAOTDDFT
from gpaw.lcaotddft.dipolemomentwriter import DipoleMomentWriter
from gpaw.lcaotddft.densitymatrix import DensityMatrix
from gpaw.lcaotddft.frequencydensitymatrix import FrequencyDensityMatrix
from gpaw.tddft.units import au_to_as
from gpaw.tddft.folding import frequencies
from gpaw.utilities.timelimit import TimeLimiter, time_to_seconds
from parallel_util import get_parallel
import numpy as np
from ase.units import Hartree, Bohr
from gpaw.external import ConstantElectricField
from gpaw.lcaotddft.laser import GaussianPulse
from gpaw.lcaotddft.wfwriter import WaveFunctionWriter
from gpaw.lcaotddft.restartfilewriter import RestartFileWriter

fxc='RPA'
# Read converged ground-state file
td_calc = LCAOTDDFT('gs.gpw', fxc=fxc, parallel=get_parallel(world))
# Attach any data recording or analysis tools
DipoleMomentWriter(td_calc, 'dm.dat')
WaveFunctionWriter(td_calc, 'wf.ulm')
# Use 'td.gpw' as restart file
RestartFileWriter(td_calc, 'td.gpw')
# Kick
td_calc.absorption_kick([1e-5, 0.0, 0.0])
# Propagate
td_calc.propagate(20, 1500)
# Save the state for restarting later
td_calc.write('td.gpw', mode='all')
