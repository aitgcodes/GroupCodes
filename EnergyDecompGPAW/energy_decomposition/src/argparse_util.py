import os
import argparse
try:
    from gpaw.mpi import world
except ImportError:
    class Communicator():
        def __init__(self, size, rank):
            self.size = size
            self.rank = rank

        def barrier(self):
            pass

    world = Communicator(size=1, rank=0)


def create_directory_path(dpath):
    if dpath == '':
        return
    if world.rank == 0:
        if not os.path.isdir(dpath):
            os.makedirs(dpath)
    world.barrier()


def DirectoryPathType(dpath):
    create_directory_path(dpath)
    return dpath


def FilePathType(fpath):
    create_directory_path(os.path.dirname(fpath))
    return fpath


def ExistingPathType(path):
    if not os.path.exists(path):
        raise ValueError(f'Given path does not exist: {path}')
    return path


def FloatOrStrType(value):
    try:
        return float(value)
    except ValueError:
        return value


def IntOrStrType(value):
    try:
        return int(value)
    except ValueError:
        return value


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paths', nargs='+', type=ExistingPathType)
    parser.add_argument('--fpaths', nargs='+', type=FilePathType)
    parser.add_argument('--dpaths', nargs='+', type=DirectoryPathType)
    parser.add_argument('--values', nargs='+', type=FloatOrStrType)
    args = parser.parse_args()

    if args.paths is not None:
        print(args.paths)
    if args.fpaths is not None:
        print(args.fpaths)
    if args.dpaths is not None:
        print(args.dpaths)
    if args.values is not None:
        for val in args.values:
            print(val, type(val))
