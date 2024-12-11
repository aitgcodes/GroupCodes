import unittest
import guesthost as gh
import numpy as np
from guesthost.analysis.functions_modules import print_file

class TestFragmentize(unittest.TestCase):

    def test_create_fragments(self):
        # Prepare data and inputs for the test
        mpb_sys_ini_path = "structures/mpb.xyz"
        trj = gh.Trajectory(mpb_sys_ini_path, order=True)
        
        ma_inds = [0, 1, 2, 3, 4, 5, 6, 7]
        ma_lists = trj.fragmentize(64, ma_inds, gh.MethylAmmonium, ordering="unitcell")
        
        host_inds = [8, 9, 10, 11]
        host_lists = trj.fragmentize(64, host_inds, gh.Host, ordering="unitcell")

        # Add assertions to verify the results
        self.assertIsInstance(ma_lists, list)
        self.assertIsInstance(host_lists, list)
        self.assertGreater(len(ma_lists), 0)
        self.assertGreater(len(host_lists), 0)

    def test_invalid_input(self):
        mpb_sys_ini_path = "path/to/invalid/file.xyz"
        with self.assertRaises(Exception):
            gh.Trajectory(mpb_sys_ini_path)

if __name__ == '__main__':
    unittest.main()

