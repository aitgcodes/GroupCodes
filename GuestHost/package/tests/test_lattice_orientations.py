import numpy as np
import unittest
import guesthost as gh
from guesthost.analysis.functions_modules import print_file

class TestLatticeOrientations(unittest.TestCase):
    def setUp(self):
        # Set up the test with required inputs
        self.mpb_sys_ini_path = "structures/mpb_trajectory.xyz"  # Path to your input file
        self.trj = gh.Trajectory(self.mpb_sys_ini_path, order=True)
        self.n = 4  # Number of unitcells in each direction
        self.all_inds = [320, 384, 385, 386, 256, 387, 388, 389, 0, 64, 65, 66]
        self.ma_inds = [0, 1, 2, 3, 4, 5, 6, 7]
        self.host_inds = [8, 9, 10, 11]

    def test_lattice_orientations(self):
        # Fragmentize and get orientation data
        ma_lists = self.trj.fragmentize(64, self.ma_inds, gh.MethylAmmonium, ordering="unitcell")
        host_lists = self.trj.fragmentize(64, self.host_inds, gh.Host, ordering="unitcell")
        lattices = self.trj.create_lattice(
            64, [self.ma_inds, self.host_inds], [gh.MethylAmmonium, gh.Host],
            (self.n, self.n, self.n), gh.HPLattice, ordering="type", unit_order=self.all_inds
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
        ref_t = np.loadtxt('reference_data/reference_thetas.dat')
        ref_p = np.loadtxt('reference_data/reference_phis.dat')

        # Compare with reference data
        np.testing.assert_allclose(t_dat, ref_t, atol=1e-6, err_msg="Theta data does not match reference.")
        np.testing.assert_allclose(p_dat, ref_p, atol=1e-6, err_msg="Phi data does not match reference.")

if __name__ == "__main__":
    unittest.main()
