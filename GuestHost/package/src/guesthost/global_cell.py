def get_global_index_data(lattice):
    return lattice.get_global_index_data()


def global_cell_vectors_OR(lattice, dir=0, global_data=None):
    return lattice.global_cell_vectors_OR(dir=dir, global_data=global_data)


def global_cell_vectors_OR_fromcell(lattice, dir=0):
    return lattice.global_cell_vectors_OR_fromcell(dir=dir)


def global_cell_lengths_angles_OR(lattice, dir=0, global_data=None, fromcell=False):
    return lattice.global_cell_lengths_angles_OR(
        dir=dir,
        global_data=global_data,
        fromcell=fromcell,
    )


def global_cell_lengths_OR(lattice, dir=0, global_data=None, fromcell=False):
    return lattice.global_cell_lengths_OR(
        dir=dir,
        global_data=global_data,
        fromcell=fromcell,
    )


def global_cell_angles_OR(lattice, dir=0, global_data=None, fromcell=False):
    return lattice.global_cell_angles_OR(
        dir=dir,
        global_data=global_data,
        fromcell=fromcell,
    )
