import os
import sys
sys.path.insert(0, os.path.dirname(__file__))


def spec(dm_fpath, spec_fpath):
    from gpaw.tddft.spectrum import photoabsorption_spectrum

    # Basic settings
    folding = 'Gauss'
    width = 0.07

    # Calculate spectrum
    photoabsorption_spectrum(dm_fpath, spec_fpath, folding, width, 0, 10, 0.01)


if __name__ == '__main__':
    import argparse
    from argparse_util import FilePathType

    parser = argparse.ArgumentParser()
    parser.add_argument('dm_fpath')
    parser.add_argument('spec_fpath', type=FilePathType)
    args = parser.parse_args()

    spec(args.dm_fpath, args.spec_fpath)
