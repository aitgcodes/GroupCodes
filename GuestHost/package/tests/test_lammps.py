from pathlib import Path

import numpy as np

import guesthost as gh


def _traj_path():
    return Path(__file__).resolve().parent / "structures" / "lammps_trajectory" / "trajectory.lmp"


def test_load_lammps_trajectory_metadata_and_compute():
    ltraj = gh.load_lammps_trajectory(
        _traj_path(),
        unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
        supercell_size=(4, 4, 4),
        steps=range(1, 3),
    )

    assert isinstance(ltraj, gh.LatticeTrajectory)
    assert len(ltraj.lattices) == 2
    assert ltraj.steps == [4260000, 4360000]
    assert ltraj.times == [426.0, 436.0]
    assert ltraj.data["stepnumber"] == ltraj.steps

    theta_phi = ltraj.compute(gh.HPLattice.ucell_theta_phi, dir=0)
    assert len(theta_phi) == 2
    assert len(theta_phi[0]) == 64
    assert len(theta_phi[0][0]) == 2

    xi_ma = ltraj.compute_layerwise(gh.global_xi_MA, dir=0)
    assert len(xi_ma) == 2
    assert all(np.isfinite(x) for x in xi_ma)


def test_load_lammps_trajectory_timestep_fallback():
    ltraj = gh.load_lammps_trajectory(
        _traj_path(),
        unitdata=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
        supercell_size=(4, 4, 4),
        steps=1,
        timestep_fs=0.2,
        prefer_logtime=False,
    )

    assert ltraj.steps == [4260000]
    assert ltraj.times == [852.0]


def test_load_lammps_trajectories_from_directory(tmp_path):
    source_dir = _traj_path().parent
    run_dir = tmp_path / "01_npt1.1gpa"
    run_dir.mkdir()
    for name in ("trajectory.lmp", "log.lammps"):
        (run_dir / name).write_bytes((source_dir / name).read_bytes())

    trajectories = gh.load_lammps_trajectories(
        tmp_path,
        unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
        supercell_size=(4, 4, 4),
        mddir_regex=r"^01",
        dump_regex=r"trajectory\.lmp$",
        steps=1,
    )

    assert list(trajectories) == ["1.1"]
    assert isinstance(trajectories["1.1"], gh.LatticeTrajectory)
