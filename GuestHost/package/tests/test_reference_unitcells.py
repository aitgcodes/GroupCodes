import json
from pathlib import Path

import numpy as np
import guesthost as gh


def _trajectory():
    traj_path = Path(__file__).parent / "structures" / "mpb_trajectory.xyz"
    return gh.Trajectory(traj_path, order=True)


def _manual_lattices(trj):
    all_inds = [320, 384, 385, 386, 256, 387, 388, 389, 0, 64, 65, 66]
    ma_inds = [0, 1, 2, 3, 4, 5, 6, 7]
    host_inds = [8, 9, 10, 11]
    return trj.create_lattice(
        64,
        [ma_inds, host_inds],
        [gh.MethylAmmonium, gh.Host],
        (4, 4, 4),
        gh.HPLattice,
        ordering="type",
        unit_order=all_inds,
    )


def _load_unitcell_json(name):
    data_path = Path(gh.__file__).resolve().parent / "data" / "unitcells" / name
    with open(data_path, "r") as f:
        data = json.load(f)

    out = []
    for entry in data:
        converted = dict(entry)
        converted["pb_axis"] = (converted["pb_axis"][0], converted["pb_axis"][1])
        out.append(converted)
    return out


def test_reference_unit_order_constant():
    assert len(gh.UNITCELL_INDEXDATA_MPB_4x4x4) == 64
    assert len(gh.UNITCELL_INDEXDATA_MPB_8x8x8) == 512
    assert len(gh.MPB_SYS_4x4x4) == 768
    assert len(gh.MPB_SYS_8x8x8) == 6144
    assert gh.UNIT_ORDER_MPB_4x4x4 == [320, 384, 385, 386, 256, 387, 388, 389, 0, 64, 65, 66]
    assert gh.UNIT_ORDERS_MPB_4x4x4[0] == gh.UNIT_ORDER_MPB_4x4x4
    assert gh.UNIT_ORDER_MPB_8x8x8 == [320, 384, 385, 386, 256, 387, 388, 389, 0, 64, 65, 66]
    assert gh.UNIT_ORDERS_MPB_8x8x8[0] == gh.UNIT_ORDER_MPB_8x8x8


def test_reference_unitcell_data_is_loaded_from_json():
    expected_keys = {
        "c_ind",
        "n_ind",
        "h_inds_c",
        "h_inds_n",
        "pb_inds",
        "br_inds",
        "pb_axis",
        "br_axis",
        "br_axis_neg",
    }
    assert isinstance(gh.UNITCELL_INDEXDATA_MPB_4x4x4[0]["pb_axis"], tuple)
    assert isinstance(gh.UNITCELL_INDEXDATA_MPB_8x8x8[0]["pb_axis"], tuple)
    assert set(gh.UNITCELL_INDEXDATA_MPB_4x4x4[0]) == expected_keys
    assert set(gh.UNITCELL_INDEXDATA_MPB_8x8x8[0]) == expected_keys
    assert gh.UNITCELL_INDEXDATA_MPB_8x8x8[-1]["c_ind"] == 5759
    assert gh.UNITCELL_INDEXDATA_MPB_8x8x8[-1]["n_ind"] == 5695


def test_unitcell_detector_matches_cached_json_4x4x4():
    expected = _load_unitcell_json("unitcell_indexdata_mpb_4x4x4.json")
    computed = gh.get_unitcell_indexdata(gh.MPB_SYS_4x4x4)
    assert computed == expected
    assert gh.UNITCELL_INDEXDATA_MPB_4x4x4 == expected


def test_unitcell_detector_can_include_ortho_axes():
    computed = gh.get_unitcell_indexdata(gh.MPB_SYS_4x4x4, ortho_axes=True)
    first = computed[0]

    assert "pb_ortho_axes" in first
    assert "pb_ortho_axes_small" in first
    assert first["pb_ortho_axes_small"] == (0, [[16, 7, 5], [4, 49, 17], [1, 20, 52]])
    assert first["pb_ortho_axes"] == (0, [[32, 7, 5], [8, 49, 17], [2, 20, 52]])


def test_reference_lattice_builder_matches_manual_lattice():
    trj = _trajectory()
    manual = _manual_lattices(trj)[1]
    inferred = trj.create_hplattice(
        unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
        supercell_size=(4, 4, 4),
    )[1]

    manual_theta, manual_phi = manual.all_ucell_theta_phi(dir=0)
    inferred_theta, inferred_phi = inferred.all_ucell_theta_phi(dir=0)
    np.testing.assert_allclose(manual_theta, inferred_theta)
    np.testing.assert_allclose(manual_phi, inferred_phi)
    np.testing.assert_allclose(manual.all_ucell_eta_MA(dir=0), inferred.all_ucell_eta_MA(dir=0))
    np.testing.assert_allclose(manual.all_ucell_eta_lat(dir_coup=0), inferred.all_ucell_eta_lat(dir_coup=0))


def test_reference_lattice_builder_matches_orientation_reference_data():
    trj = _trajectory()
    lattices = trj.create_hplattice(
        unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
        supercell_size=(4, 4, 4),
    )
    lattraj = gh.LatticeTrajectory(lattices)
    dat = lattraj.compute_for_all(gh.HPLattice.all_ucell_theta_phi, 0)

    theta = []
    phi = []
    for t, p in dat[1:]:
        theta.append(t.flatten())
        phi.append(p.flatten())
    theta = np.array(theta)
    phi = np.array(phi)

    ref_dir = Path(__file__).parent / "reference_data"
    ref_theta = np.loadtxt(ref_dir / "reference_thetas.dat")
    ref_phi = np.loadtxt(ref_dir / "reference_phis.dat")

    np.testing.assert_allclose(theta, ref_theta, atol=1e-6)
    np.testing.assert_allclose(phi, ref_phi, atol=1e-6)
