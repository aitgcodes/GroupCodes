import numpy as np
import unittest
from pathlib import Path
import guesthost as gh

class TestUnitcellFunctions(unittest.TestCase):
    def setUp(self):
        # Load the test structure, following test_lattice_orientations.py pattern
        self.mpb_sys_ini_path = Path(__file__).parent / "structures" / "mpb_trajectory.xyz"
        self.trj = gh.Trajectory(self.mpb_sys_ini_path, order=True)
        self.n = 4  # Number of unitcells in each direction

        lattices = self.trj.create_hplattice(
            unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
            supercell_size=(self.n, self.n, self.n),
        )
        self.lattraj = gh.LatticeTrajectory(lattices)

        # Get the first system and lattice for testing
        self.sys = self.trj.atoms_list[0]
        self.lattice = lattices[0]
        self.ind = (0, 0, 0)

    def test_ucell_localcellvecs(self):
        """Test local cell vectors computation"""
        vecs = self.lattice.ucell_localcellvecs(self.ind)
        self.assertEqual(len(vecs), 3)
        for v in vecs:
            self.assertEqual(len(v), 3)
            self.assertTrue(np.all(np.isfinite(v)))

    def test_ucell_localcellparameters(self):
        """Test local cell parameters computation"""
        params = self.lattice.ucell_localcellparameters(self.ind)
        required_keys = ['a', 'b', 'c', 'α', 'β', 'γ', 'vol']
        for key in required_keys:
            self.assertIn(key, params)
            self.assertTrue(np.isfinite(params[key]))

    def test_ucell_theta_phi(self):
        """Test theta phi computation"""
        theta, phi = self.lattice.ucell_theta_phi(self.ind, 0)
        self.assertTrue(np.isfinite(theta))
        self.assertTrue(np.isfinite(phi))
        self.assertGreaterEqual(theta, 0)
        self.assertLessEqual(theta, 180)

    def test_ucell_η_MA(self):
        """Test MA order parameter computation"""
        eta = self.lattice.ucell_η_MA(self.ind, 0)
        self.assertEqual(len(eta), 2)
        self.assertTrue(np.all(np.isfinite(eta)))

    def test_compute_for_all(self):
        """Test computing unitcell functions for all frames using LatticeTrajectory"""
        # Compute theta phi for all frames
        results = self.lattraj.compute_for_all(lambda lat: lat.ucell_theta_phi(self.ind, 0))
        self.assertEqual(len(results), self.trj.nt)
        for res in results:
            theta, phi = res
            self.assertTrue(np.isfinite(theta))
            self.assertTrue(np.isfinite(phi))
            self.assertGreaterEqual(theta, 0)
            self.assertLessEqual(theta, 180)

    def test_all_ucell_theta_phi(self):
        theta, phi = self.lattice.all_ucell_theta_phi(dir=0)
        self.assertEqual(theta.shape, (self.n, self.n, self.n))
        self.assertEqual(phi.shape, (self.n, self.n, self.n))
        self.assertTrue(np.all(np.isfinite(theta)))
        self.assertTrue(np.all(np.isfinite(phi)))

    def test_U_r_and_ucell_ω(self):
        """Test rotation matrix and ω computation"""
        eta = np.array([0.7, 0.3])
        U = self.lattice.U_r([0, 0], [np.pi, np.pi])
        self.assertEqual(U.shape, (2, 2))
        omega = self.lattice.ucell_ω(eta, [0, 0])
        self.assertEqual(len(omega), 2)
        self.assertTrue(np.allclose(omega, eta, atol=1e-10))  # For R=[0,0], U=I

if __name__ == "__main__":
    unittest.main()
