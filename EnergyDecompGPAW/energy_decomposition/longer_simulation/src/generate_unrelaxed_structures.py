from ase.io import write
from ase.cluster.icosahedron import Icosahedron
from ase.cluster.octahedron import Octahedron


def main():
    symbol = 'Ag'

    for n in range(3, 7):
        atoms = Icosahedron(symbol, n)
        write('unrlx-ico-{}{}.xyz'.format(symbol, len(atoms)), atoms)

    for cutoff in range(3, 6):
        length = 2 * cutoff + 1
        atoms = Octahedron(symbol, length=length, cutoff=cutoff)
        write('unrlx-cubo-{}{}.xyz'.format(symbol, len(atoms)), atoms)

    for cutoff in range(2, 4):
        length = 3 * cutoff + 1
        atoms = Octahedron(symbol, length=length, cutoff=cutoff)
        write('unrlx-rto-{}{}.xyz'.format(symbol, len(atoms)), atoms)


if __name__ == '__main__':
    main()
