from ase.io import read


def load_system_from_xyz(file_path, **kwargs):
    """Load an atomic system from an XYZ file with ASE."""
    return read(file_path, **kwargs)
