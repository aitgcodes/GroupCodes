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


def test_layerwise_lattice_order_shapes_and_consistency():
    lattice = _reference_lattice()
    cellshape = (4, 4, 4)

    omega_layers = gh.layerwise_OmegaLat(lattice, dir_coup=0, cellshape=cellshape)
    assert len(omega_layers) == cellshape[0]
    assert all(om.shape == (2,) for om in omega_layers)

    eta_layers = gh.layerwise_eta_Lat(lattice, dir_coup=0, cellshape=cellshape)
    assert len(eta_layers) == cellshape[0]
    assert all(layer.shape == (4, 4) for layer in eta_layers)

    xi_layers = gh.global_layerwise_xi(lattice, dir_coup=0, cellshape=cellshape)
    np.testing.assert_allclose(
        xi_layers,
        [0.5 * (om[0] - om[1]) for om in omega_layers],
    )

    xi_global = gh.global_xi(lattice, dir_coup=0, cellshape=cellshape)
    np.testing.assert_allclose(xi_global, np.mean(xi_layers))

    S = gh.global_S(lattice, dir_coup=0, cellshape=cellshape)
    assert len(S) == cellshape[0]
    assert S[0] >= 0.0


def test_layerwise_alloctahedra_and_ma_coupling_shapes():
    lattice = _reference_lattice()
    cellshape = (4, 4, 4)

    eta_all = gh.global_layerwise_eta_alloctahedra(
        lattice, dir_coup=1, cellshape=cellshape
    )
    xi_all = gh.global_layerwise_xi_alloctahedra(
        lattice, dir_coup=1, cellshape=cellshape
    )
    xi_new = gh.global_layerwise_xi_new(lattice, dir_coup=1, cellshape=cellshape)
    assert len(eta_all) == cellshape[1]
    assert len(xi_all) == cellshape[1]
    assert len(xi_new) == cellshape[1]

    omega_ma = gh.layerwise_OmegaMA(lattice, dir=2, cellshape=cellshape)
    eta_ma = gh.layerwise_eta_MA(lattice, dir=2, cellshape=cellshape)
    xi_ma = gh.global_xi_MA(lattice, dir=2, cellshape=cellshape)
    omega_c = gh.layerwise_OmegaCoupling(lattice, dir=2, cellshape=cellshape)
    xi_c = gh.global_xi_coupling(lattice, dir=2, cellshape=cellshape)

    assert len(omega_ma) == cellshape[2]
    assert len(eta_ma) == cellshape[2]
    assert isinstance(xi_ma, float)
    assert len(omega_c) == cellshape[2]
    assert isinstance(xi_c, float)
    assert np.isfinite(gh.global_volume_local(lattice))


def test_compute_all_ucells_helpers():
    lattice = _reference_lattice()
    cellshape = (4, 4, 4)

    vals = gh.compute_all_ucells(
        lambda lat, ind: lat.ucell_udata(ind)["c_ind"],
        lattice,
    )
    assert len(vals) == 64

    mat = gh.compute_all_ucells_matrix(
        lambda lat, ind: np.array(lat.ucell_localcelllengths(ind)),
        lattice,
        cellshape=cellshape,
    )
    assert mat.shape == cellshape
    assert mat[0, 0, 0].shape == (3,)

    data = gh.compute_all_ucells_data(
        lambda lat, ind: lat.ucell_theta_phi(ind, dir=0),
        lattice,
        cellshape=cellshape,
    )
    assert data.shape == cellshape
    assert len(data[0, 0, 0]) == 2
