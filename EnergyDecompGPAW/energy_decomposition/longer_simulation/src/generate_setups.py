from gpaw.atom.generator import Generator
from gpaw.atom.basis import BasisMaker
from gpaw.atom.configurations import parameters_extra


def main():
    atom = 'Ag'
    name = 'my'
    args = parameters_extra[atom]  # Choose the 11-electron setup
    args.update(dict(name=name, use_restart_file=False, exx=True))

    # Setup and basis for relaxation (PBE)
    xc = 'PBE'

    # Generate setup
    generator = Generator(atom, xc, scalarrel=True)
    generator.run(write_xml=True, **args)

    # Generate basis
    bm = BasisMaker(atom, name='{}.{}'.format(name, xc), xc=xc, run=False)
    bm.generator.run(write_xml=False, **args)
    basis = bm.generate(zetacount=2, polarizationcount=1,
                        jvalues=[0, 1],
                        )
    basis.write_xml()

    # Setup and basis for TDDFT (GLLBSC)
    xc = 'GLLBSC'

    # Generate setup
    generator = Generator(atom, xc, scalarrel=True)
    generator.run(write_xml=True, **args)

    # Generate basis
    bm = BasisMaker(atom, name='{}.{}.p'.format(name, xc), xc=xc, run=False)
    bm.generator.run(write_xml=False, **args)
    basis = bm.generate(zetacount=2, polarizationcount=0,
                        jvalues=[0, 1, 2],
                        )
    basis.write_xml()


if __name__ == '__main__':
    main()
