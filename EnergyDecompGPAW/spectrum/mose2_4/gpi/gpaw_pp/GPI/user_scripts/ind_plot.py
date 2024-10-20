# A sample script to generate plots of induced quantities 
# Adapted from scripts on GPAW website
import numpy as np
import matplotlib.pyplot as plt

from ase.io import read

syst = "Na8"

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
    plt.figure(figsize=(8, 3.5))
    ax = plt.subplot(1, 1, 1)
    X, Y = np.meshgrid(x, y)
    dmax = max(d_yx.min(), d_yx.max())
    vmax = 0.9 * dmax
    vmin = -vmax
    plt.pcolormesh(X, Y, d_yx.T, cmap='RdBu_r', vmin=vmin, vmax=vmax,shading='auto')
    contours = np.sort(np.outer([-1, 1], [0.02]).ravel() * dmax)
    plt.contour(X, Y, d_yx.T, contours, cmap='RdBu_r', vmin=-1e-10, vmax=1e-10)
    for atom in atoms:
        pos = atom.position
        plt.scatter(pos[0], pos[1], s=100, c='k', marker='o')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim([x[0], x[-1]])
    plt.ylim([y[0], y[-1]])
    ax.set_aspect('equal')

    plt.title(f'{property} of {syst} at {freq:.2f} eV')
    plt.tight_layout()
    plt.savefig(f'{fstr}_{freq:.2f}.png')

do(1.12,'den')
do(1.12,'pot')
do(1.12,'feh')
do(2.48,'den')
do(2.48,'feh')
