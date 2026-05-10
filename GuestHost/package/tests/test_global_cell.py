from importlib.resources import files

import numpy as np

import guesthost as gh


def _reference_lattice():
    data_path = files("guesthost").joinpath("data", "structures", "mpb_cubic_4x4x4.xyz")
    trj = gh.Trajectory(str(data_path), order=True)
    return trj.create_hplattice(
        supercell_size=(4, 4, 4),
        unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
    )[0]


def test_global_index_data_uses_pb_indices():
    lattice = _reference_lattice()
    global_data = lattice.get_global_index_data()

    assert global_data["supercell_size"] == 4
    assert len(global_data["axes_OR"]) == 3

    pb_indices = {
        lattice.get_motif(ind).host.indices[0]
        for ind in np.ndindex(lattice.nx, lattice.ny, lattice.nz)
    }
    for origin, endpt in global_data["axes_OR"]:
        assert origin in pb_indices
        assert endpt in pb_indices


def test_global_cell_vectors_lengths_and_angles():
    lattice = _reference_lattice()

    for direction in range(3):
        vecs = gh.global_cell_vectors_OR(lattice, dir=direction)
        assert len(vecs) == 3
        assert all(np.asarray(v).shape == (3,) for v in vecs)

        lengths, angles = gh.global_cell_lengths_angles_OR(lattice, dir=direction)
        np.testing.assert_allclose(lengths, [np.linalg.norm(v) for v in vecs])
        assert len(angles) == 3
        assert all(0.0 <= angle <= 180.0 for angle in angles)

        np.testing.assert_allclose(
            gh.global_cell_lengths_OR(lattice, dir=direction),
            lengths,
        )
        np.testing.assert_allclose(
            gh.global_cell_angles_OR(lattice, dir=direction),
            angles,
        )


def test_global_cell_vectors_fromcell_definition():
    lattice = _reference_lattice()
    cell = np.asarray(lattice.cell)

    vecs = lattice.global_cell_vectors_OR_fromcell(dir=0)
    np.testing.assert_allclose(vecs[0], 2.0 * cell[0])
    np.testing.assert_allclose(vecs[1], cell[1] - cell[2])
    np.testing.assert_allclose(vecs[2], cell[1] + cell[2])

    lengths, angles = lattice.global_cell_lengths_angles_OR(dir=0, fromcell=True)
    np.testing.assert_allclose(lengths, [np.linalg.norm(v) for v in vecs])
    np.testing.assert_allclose(angles, [90.0, 90.0, 90.0], atol=1e-12)
