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


# Temporal shape of the time-dependent potential
pulse = GaussianPulse(1e-5, 10e3, 4.14, 0.3, 'sin')
# Spatial shape of the time-dependent potential
ext = ConstantElectricField(Hartree / Bohr, [0, 0, 1])
# Full time-dependent potential
td_potential = {'ext': ext, 'laser': pulse}

# Write the temporal shape to a file
pulse.write('pulse.dat', np.arange(0, 30e3, 5.0))

fxc='RPA'
# Set up the time-propagation calculation
td_calc = LCAOTDDFT('gs.gpw', fxc=fxc, parallel=get_parallel(world),
                     td_potential=td_potential,
                     txt='tdpulse.out')

# Attach the data recording and analysis tools
DipoleMomentWriter(td_calc, 'dmpulse.dat')
#WaveFunctionWriter(td_calc, 'wf.ulm')

# Use 'td.gpw' as restart file
RestartFileWriter(td_calc, 'td.gpw')

# Propagate
td_calc.propagate(20, 6000)

# Save the state for restarting later
td_calc.write('td.gpw', mode='all')
