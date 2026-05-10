from pathlib import Path

import numpy as np
import guesthost as gh


def _resolve_test_paths():
    pkg_root = Path(gh.__file__).resolve().parents[2]
    tests_root = pkg_root / "tests"
    data_path = tests_root / "structures" / "mpb_trajectory.xyz"
    ref_dir = tests_root / "reference_data"
    if not data_path.exists():
        raise FileNotFoundError(f"Missing test structure at {data_path}")
    if not ref_dir.exists():
        raise FileNotFoundError(f"Missing reference data at {ref_dir}")
    return data_path, ref_dir


def main():
    data_path, ref_dir = _resolve_test_paths()

    trj = gh.Trajectory(str(data_path), order=True)
    lattices = trj.create_hplattice(
        unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
        supercell_size=(4, 4, 4),
    )
    lattraj = gh.LatticeTrajectory(lattices)

    ax = 0
    dat = lattraj.compute_for_all(gh.HPLattice.all_ma_orientations, ax)

    t_dat = []
    p_dat = []
    for val in dat[1:]:
        t, p = val
        t_dat.append(t.flatten())
        p_dat.append(p.flatten())
    t_dat = np.array(t_dat)
    p_dat = np.array(p_dat)

    ref_t = np.loadtxt(ref_dir / "reference_thetas.dat")
    ref_p = np.loadtxt(ref_dir / "reference_phis.dat")

    np.testing.assert_allclose(t_dat, ref_t, atol=1e-6, err_msg="Theta data does not match reference.")
    np.testing.assert_allclose(p_dat, ref_p, atol=1e-6, err_msg="Phi data does not match reference.")
    print("guesthost lattice orientation check: OK")


if __name__ == "__main__":
    main()
