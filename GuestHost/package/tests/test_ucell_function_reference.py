import json
from pathlib import Path

import numpy as np
import pytest

import guesthost as gh


DIRS_PERP = [(1, 2), (2, 0), (0, 1)]


@pytest.fixture(scope="module")
def reference_lattice():
    structure = Path(__file__).parent / "structures" / "1.1gpa_step10_after150ps_allatoms.xyz"
    trj = gh.Trajectory(structure, order=True)
    return trj.create_hplattice(
        unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
        supercell_size=(4, 4, 4),
    )[0]


@pytest.fixture(scope="module")
def reference_data():
    ref_path = Path(__file__).parent / "reference_data" / "julia_reference_ucell_functions.json"
    with open(ref_path, "r") as f:
        return json.load(f)


def _assert_close(actual, expected, atol=1e-8):
    if isinstance(expected, dict):
        assert set(actual) == set(expected)
        for key in expected:
            _assert_close(actual[key], expected[key], atol=atol)
        return

    if isinstance(expected, list) and expected and isinstance(expected[0], dict):
        assert len(actual) == len(expected)
        for actual_item, expected_item in zip(actual, expected):
            _assert_close(actual_item, expected_item, atol=atol)
        return

    np.testing.assert_allclose(actual, expected, atol=atol)


def _motif_index(flat_index):
    return np.unravel_index(flat_index, (4, 4, 4))


def _localcellparameters_for_reference(lattice, ind):
    params = lattice.ucell_localcellparameters(ind)
    return {
        "a": params["a"],
        "b": params["b"],
        "c": params["c"],
        "α": params["alpha"],
        "β": params["beta"],
        "γ": params["gamma"],
        "vol": params["vol"],
    }


def _hostrelative_values(lattice, ind, method_name):
    method = getattr(lattice, method_name)
    return {
        str(dir_coup + 1): [
            method(ind, dir=dir, dir_coup=dir_coup)
            for dir in DIRS_PERP[dir_coup]
        ]
        for dir_coup in range(3)
    }


def _plane_distance_values(lattice, ind):
    return {
        str(dir_coup + 1): [
            lattice.ucell_br_distance_from_plane(ind, dir_pl=[dir, dir_coup], dir_br=dir)
            for dir in DIRS_PERP[dir_coup]
        ]
        for dir_coup in range(3)
    }


@pytest.mark.parametrize("flat_index", [0, 21, 63])
def test_ucell_geometry_matches_reference(reference_lattice, reference_data, flat_index):
    ind = _motif_index(flat_index)
    ref = reference_data[flat_index]

    _assert_close(reference_lattice.ucell_localcellvecs(ind), ref["localcellvecs"])
    _assert_close(_localcellparameters_for_reference(reference_lattice, ind), ref["localcellparameters"])
    _assert_close(reference_lattice.ucell_localcelllengths(ind), ref["localcelllengths"])
    _assert_close(reference_lattice.ucell_localcellangles(ind), ref["localcellangles"])
    _assert_close(reference_lattice.ucell_localcell_lengths_angles(ind), ref["localcelllengths_angles"])

    localcellvecs_ortho = {
        f"dir_{dir + 1}": reference_lattice.ucell_localcellvecs_ortho(ind, dir=dir)
        for dir in range(3)
    }
    localcellvecs_ortho_small = {
        f"dir_{dir + 1}": reference_lattice.ucell_localcellvecs_ortho(ind, dir=dir, small=True)
        for dir in range(3)
    }
    localcelllengths_ortho = {
        f"dir_{dir + 1}": reference_lattice.ucell_localcelllengths_ortho(ind, dir=dir)
        for dir in range(3)
    }
    localcellangles_ortho = {
        f"dir_{dir + 1}": reference_lattice.ucell_localcellangles_ortho(ind, dir=dir)
        for dir in range(3)
    }
    localcelllengths_angles_ortho = {
        f"dir_{dir + 1}": reference_lattice.ucell_localcell_lengths_angles_ortho(ind, dir=dir)
        for dir in range(3)
    }

    _assert_close(localcellvecs_ortho, ref["localcellvecs_ortho"])
    _assert_close(localcellvecs_ortho_small, ref["localcellvecs_ortho_small"])
    _assert_close(localcelllengths_ortho, ref["localcelllengths_ortho"])
    _assert_close(localcellangles_ortho, ref["localcellangles_ortho"])
    _assert_close(localcelllengths_angles_ortho, ref["localcelllengths_angles_ortho"])


@pytest.mark.parametrize("flat_index", [0, 21, 63])
def test_ucell_ma_orientation_matches_reference(reference_lattice, reference_data, flat_index):
    ind = _motif_index(flat_index)
    ref = reference_data[flat_index]

    theta_phi = [reference_lattice.ucell_theta_phi_dict(ind, dir=dir) for dir in range(3)]
    theta_phi_withcutoff = [
        reference_lattice.ucell_theta_phi_withcutoff(ind, dir=dir)
        for dir in range(3)
    ]
    xyz_from_theta_phi = [
        reference_lattice.ucell_xyz(tp["theta"], tp["phi"], theta_dir=dir)
        for dir, tp in enumerate(theta_phi)
    ]

    _assert_close(theta_phi, ref["theta_phi"])
    _assert_close(theta_phi_withcutoff, ref["theta_phi_withcutoff"])
    _assert_close(xyz_from_theta_phi, ref["xyz_from_theta_phi"])
    _assert_close([reference_lattice.ucell_eta_MA(ind, dir=dir) for dir in range(3)], ref["eta_MA"])
    _assert_close([reference_lattice.ucell_eta_3d_MA(ind, dir=dir) for dir in range(3)], ref["eta_3d_MA"])


@pytest.mark.parametrize("flat_index", [0, 21, 63])
def test_ucell_host_distortions_matches_reference(reference_lattice, reference_data, flat_index):
    ind = _motif_index(flat_index)
    ref = reference_data[flat_index]

    _assert_close([reference_lattice.ucell_br_br_vectors(ind, dir=dir) for dir in range(3)], ref["br_br_vectors"])
    _assert_close([reference_lattice.ucell_br_br_distances(ind, dir=dir) for dir in range(3)], ref["br_br_distances"])
    _assert_close(
        [reference_lattice.ucell_tetragonal_lengthening_distortion_br(ind, dir=dir) for dir in range(3)],
        ref["tetragonal_lengthening_distortion_br"],
    )
    _assert_close(
        [reference_lattice.ucell_scissoring_distortion_br(ind, dir=dir) for dir in range(3)],
        ref["scissoring_distortion_br"],
    )
    _assert_close(
        [reference_lattice.ucell_scissoring_distortion_angle_br(ind, dir=dir) for dir in range(3)],
        ref["scissoring_distortion_angle_br"],
    )
    _assert_close([reference_lattice.ucell_pb_br_pb_angle(ind, dir=dir) for dir in range(3)], ref["pb_br_pb_angle"])
    _assert_close(
        _hostrelative_values(reference_lattice, ind, "ucell_pb_br_pb_angle_hostrelative"),
        ref["pb_br_pb_angle_hostrelative"],
    )
    _assert_close(_plane_distance_values(reference_lattice, ind), ref["br_distance_from_plane"])
    _assert_close(
        _hostrelative_values(reference_lattice, ind, "ucell_pb_br_pb_angle_alongaxis"),
        ref["pb_br_pb_angle_alongaxis"],
    )
    _assert_close(
        _hostrelative_values(reference_lattice, ind, "ucell_pb_br_pb_angle_alongaxis_hostrelative"),
        ref["pb_br_pb_angle_alongaxis_hostrelative"],
    )
    _assert_close([reference_lattice.ucell_eta_lat(ind, dir_coup=dir) for dir in range(3)], ref["eta_lat"])
    _assert_close(
        [reference_lattice.ucell_eta_lat_alongaxis(ind, dir_coup=dir) for dir in range(3)],
        ref["eta_lat_alongaxis"],
    )


@pytest.mark.parametrize("flat_index", [0, 21, 63])
def test_ucell_ma_and_index_helpers_match_reference(reference_lattice, reference_data, flat_index):
    ind = _motif_index(flat_index)
    ref = reference_data[flat_index]

    _assert_close(reference_lattice.ucell_MA_dihedrals(ind), ref["MA_dihedrals"])
    _assert_close(reference_lattice.ucell_ma_displacement_vecs_br(ind), ref["ma_displacement_vecs_br"])
    _assert_close(reference_lattice.ucell_CN_distance(ind), ref["CN_distance"])
    _assert_close(reference_lattice.ucell_CH_distances(ind), ref["CH_distances"])
    _assert_close(reference_lattice.ucell_NH_distances(ind), ref["NH_distances"])
    assert list(reference_lattice.ucell_pb_axis(ind)) == ref["pb_axis"]
    assert list(reference_lattice.ucell_pb_axis(ind, python_index=True)) == ref["pb_axis_python"]
    _assert_close(reference_lattice.ucell_all_octahedra_indices(ind), ref["all_octahedra_indices"])
