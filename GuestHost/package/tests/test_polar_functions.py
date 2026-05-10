from importlib.resources import files

import numpy as np

import guesthost as gh
from guesthost.polar import _polar_op_from_phi


def _reference_lattice():
    data_path = files("guesthost").joinpath("data", "structures", "mpb_cubic_4x4x4.xyz")
    trj = gh.Trajectory(str(data_path), order=True)
    return trj.create_hplattice(
        supercell_size=(4, 4, 4),
        unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
    )[0]


def test_polar_order_parameter_shapes_and_definition():
    lattice = _reference_lattice()
    cellshape = (4, 4, 4)

    result = gh.polar_order_parameter(lattice, dir_coup=0, cellshape=cellshape)
    expected_keys = {
        "dphi",
        "z",
        "z_pi",
        "z_pi_by2",
        "z_chain",
        "z_chain_pi",
        "z_chain_pi_by2",
        "S_chain",
        "S_chain_pi",
        "S_chain_pi_by2",
    }
    assert set(result) == expected_keys
    assert result["dphi"].shape == cellshape
    assert result["z"].shape == cellshape
    assert result["z_chain"].shape == (cellshape[0],)
    assert result["S_chain"].shape == (cellshape[0],)

    phi_arr = np.empty(cellshape, dtype=float)
    for ind in np.ndindex(cellshape):
        phi_arr[ind] = lattice.ucell_theta_phi(ind, dir=0)[1]
    expected = _polar_op_from_phi(phi_arr, dir_coup=0)

    for key in expected_keys:
        np.testing.assert_allclose(result[key], expected[key])


def test_polar_order_parameter_wraps_angle_differences():
    phi_arr = np.zeros((2, 1, 1))
    phi_arr[0, 0, 0] = 170.0
    phi_arr[1, 0, 0] = -170.0

    result = _polar_op_from_phi(phi_arr, dir_coup=0)
    np.testing.assert_allclose(result["dphi"][:, 0, 0], [20.0, -20.0])
    np.testing.assert_allclose(np.abs(result["z"]), 1.0)
