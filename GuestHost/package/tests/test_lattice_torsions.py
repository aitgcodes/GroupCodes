import numpy as np
import unittest
from pathlib import Path
import guesthost as gh
from guesthost.analysis.functions_modules import print_file

class TestLatticeTorsions(unittest.TestCase):
    def setUp(self):
        # Set up the test with required inputs
        self.mpb_sys_ini_path = Path(__file__).parent / "structures" / "mpb_trajectory.xyz"
        self.pref = "test"  # Prefix for file names
        self.n = 4  # Number of unitcells in each direction
        self.trj = gh.Trajectory(self.mpb_sys_ini_path)

    def test_lattice_torsions(self):
        lattices = self.trj.create_hplattice(
            unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
            supercell_size=(self.n, self.n, self.n),
        )
        lattraj = gh.LatticeTrajectory(lattices)
        
        fn = gh.HPLattice.all_ma_torsions
        ini_lat = lattices[0]
        dat_N = lattraj.compute_for_all(fn, ini_lat, htyp="N")
        dat_C = lattraj.compute_for_all(fn, ini_lat, htyp="C")

        # Collect torsion data
        t_dat_N = []
        t_dat_C = []
        for val_N, val_C in zip(dat_N, dat_C):
            t_dat_N.append(val_N.flatten())
            t_dat_C.append(val_C.flatten())
        t_dat_N = np.array(t_dat_N)
        t_dat_C = np.array(t_dat_C)

        # Load reference data
        ref_data_dir = Path(__file__).parent / "reference_data"
        ref_t_N = np.loadtxt(ref_data_dir / 'reference_Torsion_N.dat')
        ref_t_C = np.loadtxt(ref_data_dir / 'reference_Torsion_C.dat')

        # Compare with reference data
        np.testing.assert_allclose(t_dat_N, ref_t_N, atol=1e-6, err_msg="Torsion data for N does not match reference.")
        np.testing.assert_allclose(t_dat_C, ref_t_C, atol=1e-6, err_msg="Torsion data for C does not match reference.")

if __name__ == "__main__":
    unittest.main()
