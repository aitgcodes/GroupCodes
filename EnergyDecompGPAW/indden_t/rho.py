import numpy as np
from ase.io import read, write
from gpaw.lcaotddft import LCAOTDDFT
from gpaw.lcaotddft.densitymatrix import DensityMatrix
from gpaw.lcaotddft.timedensitymatrix import TimeDensityMatrix
from gpaw.lcaotddft.dipolemomentwriter import DipoleMomentWriter
from gpaw.lcaotddft.linedensity import LineDensityWriter
from gpaw.lcaotddft.linedensity import LineDensityReader
from gpaw.lcaotddft.ksdecomposition import KohnShamDecomposition

from gpaw import GPAW


# Read the ground-state file
td_calc = LCAOTDDFT('gs.gpw')
dmat = DensityMatrix(td_calc)
ksd = KohnShamDecomposition(td_calc, 'ksd.ulm')

interval=10
# Attach analysis tools
tdm=TimeDensityMatrix(td_calc, dmat, ksd, interval)

##Replay the propogation
td_calc.replay(name = 'wf.ulm', update = 'all')
