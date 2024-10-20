# A sample script to generate plots of induced quantities 
# Adapted from scripts on GPAW website
import numpy as np
import matplotlib.pyplot as plt

import sys
from ase.io import read

#syst=sys.argv[1]
freq=float(sys.argv[1])
ftyp=sys.argv[2]
print(ftyp)

def do(freq,ftyp):
    # Read cube file
    if ftyp == 'den':
       fstr = "ind"
       cube = read(f'ind_{freq:.2f}.cube', full_output=True)
       property = "Induced density"
    elif ftyp == 'pot':
       fstr = "indpot"
       cube = read(f'indpot_{freq:.2f}.cube', full_output=True)
       property = "Induced potential"
    elif ftyp == 'feh':
       fstr = "field"
       cube = read(f'field_{freq:.2f}.cube', full_output=True)
       property = "Field enhancement"
    

    #cube = read(f'{fstr}_{freq:.2f}.cube', full_output=True)

    d_g = cube['data']
    atoms = cube['atoms']
    box = np.diag(atoms.get_cell())
    ng = d_g.shape

    # Take slice of data array
    d_yx = d_g[:, :, ng[2] // 2]
    x = np.linspace(0, box[0], ng[0])
    xlabel = u'x (A)'
    y = np.linspace(0, box[1], ng[1])
    ylabel = u'y (A)'

    # Plot
    #plt.figure(figsize=(8, 8))
    #ax = plt.subplot(1, 1, 1)
    X, Y = np.meshgrid(x, y)
    #dmax = max(d_yx.min(), d_yx.max())
    #vmax = 0.9 * dmax
    #vmin = -vmax
    (dmax, dmin) = (d_yx.max(),d_yx.min())
    vmax = dmax
    vmin = dmin
    #plt.pcolormesh(X, Y, d_yx.T, cmap='RdBu_r', vmin=vmin, vmax=vmax,shading='auto')
    #contours = np.sort(np.outer([-1, 1], [0.02]).ravel() * dmax)
    #plt.contour(X, Y, d_yx.T, contours, cmap='RdBu_r', vmin=-1e-10, vmax=1e-10)
    #plt.contour(X, Y, d_yx.T, 20, cmap='RdGy')
    plt.imshow(d_yx.T, extent=[X.min(), X.max(), Y.min(), Y.max()], origin='lower', cmap='RdBu_r')
    plt.colorbar()

    sz={"Mo":180, "Se":90}
    col={"Mo":'#C54B8C', "Se":'#008000'}
    for atom in atoms:
        pos = atom.position
        plt.scatter(pos[0], pos[1], s=sz[atom.symbol], c=col[atom.symbol], marker='o')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #plt.xlim([x[0], x[-1]])
    #plt.ylim([y[0], y[-1]])
    #ax.set_aspect('equal')

    #plt.title(f'{property} of {syst} at {freq:.2f} eV')
    plt.tight_layout()
    plt.savefig(f'{fstr}_{freq:.2f}.png')

do(freq,ftyp)
