import numpy as np
import unittest
import guesthost as gh
from guesthost.analysis.functions_modules import print_file

class TestLatticeTorsions(unittest.TestCase):
    def setUp(self):
        # Set up the test with required inputs
        self.mpb_sys_ini_path = "structures/mpb_trajectory.xyz"  # Path to your input file
        self.pref = "test"  # Prefix for file names
        self.n = 4  # Number of unitcells in each direction
        self.trj = gh.Trajectory(self.mpb_sys_ini_path)
        self.all_inds = [320, 384, 385, 386, 256, 387, 388, 389, 0, 64, 65, 66]
        self.ma_inds = [0, 1, 2, 3, 4, 5, 6, 7]
        self.host_inds = [8, 9, 10, 11]

    def test_lattice_torsions(self):
        # Fragmentize and get torsion data
        lattices = self.trj.create_lattice(
            64, [self.ma_inds, self.host_inds], [gh.MethylAmmonium, gh.Host],
            (self.n, self.n, self.n), gh.HPLattice, ordering="type", unit_order=self.all_inds
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
        ref_t_N = np.loadtxt('reference_data/reference_Torsion_N.dat')
        ref_t_C = np.loadtxt('reference_data/reference_Torsion_C.dat')

        # Compare with reference data
        np.testing.assert_allclose(t_dat_N, ref_t_N, atol=1e-6, err_msg="Torsion data for N does not match reference.")
        np.testing.assert_allclose(t_dat_C, ref_t_C, atol=1e-6, err_msg="Torsion data for C does not match reference.")

if __name__ == "__main__":
    unittest.main()
