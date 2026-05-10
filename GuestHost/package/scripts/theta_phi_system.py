from pathlib import Path
import guesthost as gh

pkg_root = Path(gh.__file__).resolve().parents[2]
data_path = pkg_root / "tests" / "structures" / "mpb_trajectory.xyz"
if not data_path.exists():
    raise FileNotFoundError(f"Missing test data at {data_path}")

trj = gh.Trajectory(str(data_path), order=True)

lattices = trj.create_hplattice(
    unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
    supercell_size=(4, 4, 4),
)

lattice = lattices[1]

theta, phi = lattice.ucell_theta_phi((0, 0, 0), dir=0)
print(theta, phi)
