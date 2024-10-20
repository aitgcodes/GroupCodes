from gpaw import GPAW
from gpaw.lcaotddft.ksdecomposition import KohnShamDecomposition

# Calculate ground state with full unoccupied space
calc = GPAW('gs.gpw', txt=None).fixed_density(nbands='nao', txt='unocc.out')
calc.write('unocc.gpw', mode='all')

# Construct KS electron-hole basis
ksd = KohnShamDecomposition(calc)
ksd.initialize(calc)
ksd.write('ksd.ulm')
