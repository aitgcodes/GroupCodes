from pathlib import Path
import numpy as np
import unittest
import guesthost as gh
from guesthost.analysis.functions_modules import print_file

class TestLatticeOrientations(unittest.TestCase):
    def setUp(self):
        # Set up the test with required inputs
        self.mpb_sys_ini_path = Path(__file__).parent / "structures" / "mpb_trajectory.xyz"  # Path to your input file
        self.trj = gh.Trajectory(self.mpb_sys_ini_path, order=True)
        self.n = 4  # Number of unitcells in each direction

    def test_lattice_orientations(self):
        lattices = self.trj.create_hplattice(
            unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
            supercell_size=(self.n, self.n, self.n),
        )
        lattraj = gh.LatticeTrajectory(lattices)
        
        ax = 0  # Axis to compute orientation for
        dat = lattraj.compute_for_all(gh.HPLattice.all_ma_orientations, ax)

        # Extract theta and phi
        t_dat = []
        p_dat = []
        for val in dat[1:]:
            t, p = val
            t_dat.append(t.flatten())
            p_dat.append(p.flatten())
        t_dat = np.array(t_dat)
        p_dat = np.array(p_dat)

        # Load reference data
        ref_data_dir = Path(__file__).parent / "reference_data"
        ref_t = np.loadtxt(ref_data_dir / "reference_thetas.dat")
        ref_p = np.loadtxt(ref_data_dir / "reference_phis.dat")

        # Compare with reference data
        np.testing.assert_allclose(t_dat, ref_t, atol=1e-6, err_msg="Theta data does not match reference.")
        np.testing.assert_allclose(p_dat, ref_p, atol=1e-6, err_msg="Phi data does not match reference.")

if __name__ == "__main__":
    unittest.main()
