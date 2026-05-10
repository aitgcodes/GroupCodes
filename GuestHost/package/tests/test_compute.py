from importlib.resources import files

import numpy as np

import guesthost as gh


def _reference_lattices():
    data_path = files("guesthost").joinpath("data", "structures", "mpb_trajectory.xyz")
    trj = gh.Trajectory(str(data_path), order=True)
    return trj.create_hplattice(
        supercell_size=(4, 4, 4),
        unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
    )


def test_compute_unitcell_function_on_lattice():
    lattice = _reference_lattices()[0]

    lengths = gh.compute(lattice, gh.HPLattice.ucell_localcelllengths)
    assert len(lengths) == 64
    np.testing.assert_allclose(lengths[0], lattice.ucell_localcelllengths((0, 0, 0)))

    theta_phi = gh.compute(lattice, gh.HPLattice.ucell_theta_phi, dir=0)
    assert len(theta_phi) == 64
    assert len(theta_phi[0]) == 2


def test_compute_trajectory_and_bundle():
    lattices = _reference_lattices()
    ltraj = gh.LatticeTrajectory(lattices)

    lengths = gh.compute_trajectory(ltraj, gh.HPLattice.ucell_localcelllengths)
    assert len(lengths) == len(lattices)
    assert len(lengths[0]) == 64

    bundle = gh.compute_functions(
        ltraj,
        [
            (gh.HPLattice.ucell_theta_phi, {"dir": 0}, False),
            (gh.global_xi_MA, {"dir": 0}, True),
        ],
    )
    assert "ucell_theta_phi(dir=0)" in bundle
    assert "layerwise:global_xi_MA(dir=0)" in bundle
    assert len(bundle["ucell_theta_phi(dir=0)"]) == len(lattices)
    assert len(bundle["layerwise:global_xi_MA(dir=0)"]) == len(lattices)

    np.testing.assert_allclose(
        ltraj.compute_layerwise(gh.global_xi_MA, dir=0),
        bundle["layerwise:global_xi_MA(dir=0)"],
    )


def test_compute_default_functions_basic_keys():
    lattice = _reference_lattices()[0]
    bundle = gh.compute_default_functions(lattice, cellshape=(4, 4, 4))

    assert "ucell_theta_phi(dir=0)" in bundle
    assert "ucell_localcellparameters" in bundle
    assert "layerwise:layerwise_eta_Lat(cellshape=(4, 4, 4),dir_coup=0)" in bundle
    assert "layerwise:global_xi_coupling(cellshape=(4, 4, 4),dir=2)" in bundle
    assert len(bundle["ucell_theta_phi(dir=0)"]) == 64


def test_save_and_load_results(tmp_path):
    lattices = _reference_lattices()
    ltraj = gh.LatticeTrajectory(lattices)
    ltraj.times = [0.0, 1.0]
    ltraj.steps = [10, 20]

    results = {
        "lengths": gh.compute_trajectory(ltraj, gh.HPLattice.ucell_localcelllengths),
        "xi_ma": ltraj.compute_layerwise(gh.global_xi_MA, dir=0),
        "metadata": {"scalar": 2.0, "vector": np.array([1.0, 2.0])},
    }

    out = tmp_path / "results.h5"
    gh.save_results(out, results, trajectory=ltraj)
    loaded = gh.load_results(out)

    assert "lengths" in loaded
    assert "xi_ma" in loaded
    assert "metadata" in loaded
    np.testing.assert_allclose(loaded["time"], [0.0, 1.0])
    np.testing.assert_array_equal(loaded["steps"], [10, 20])
    np.testing.assert_allclose(loaded["xi_ma"], results["xi_ma"])
    np.testing.assert_allclose(loaded["metadata"]["vector"], [1.0, 2.0])
