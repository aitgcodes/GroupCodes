from ase import io
from ase.io.bader import attach_charges

hirsh = io.read('Hirshfeld.traj')

bader = hirsh.copy()
attach_charges(bader)

print('atom Hirshfeld Bader')
for ah, ab in zip(hirsh, bader):
    assert ah.symbol == ab.symbol
    print(f'{ah.symbol:4s} {ah.charge:9.2f} {ab.charge:5.2f}')
