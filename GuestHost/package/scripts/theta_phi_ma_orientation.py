from pathlib import Path

import guesthost as gh


def _resolve_test_path() -> Path:
    pkg_root = Path(gh.__file__).resolve().parents[2]
    data_path = pkg_root / "tests" / "structures" / "mpb_trajectory.xyz"
    if not data_path.exists():
        raise FileNotFoundError(f"Missing test data at {data_path}")
    return data_path


def main() -> None:
    data_path = _resolve_test_path()

    trj = gh.Trajectory(str(data_path), order=True)
    lattices = trj.create_hplattice(
        unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
        supercell_size=(4, 4, 4),
    )

    frame_index = 1
    lattice = lattices[frame_index]

    ax = 0
    unitcell_index = (0, 0, 0)
    theta, phi = lattice.ma_orientation(unitcell_index, ax)
    print(theta, phi)


if __name__ == "__main__":
    main()
