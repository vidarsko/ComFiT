import unittest
import sys
import os

# Adjust the path to import the comfit package
sys.path.append(os.path.abspath('../'))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

class TestBoseEinsteinCondensate(unittest.TestCase):
    def test_init_with_dimension(self):
        """Test BoseEinsteinCondensate initialization with a dimension parameter."""
        for dim in [1, 2, 3]:
            try:
                bec = cf.BoseEinsteinCondensate(dim)
                self.assertIsInstance(bec, cf.BoseEinsteinCondensate)
            except Exception as e:
                self.fail(f"Initialization failed with dimension {dim}: {e}")

    def test_dGPE_relaxer(self):
        """Test the dGPE relaxer."""
        bec = cf.BoseEinsteinCondensate(2,xRes=13,yRes=4)
        bec.conf_initial_condition_disordered()
        bec.evolve_relax(300)

        # Set the tolerance for approximation
        tolerance = 0.01

        # Check if all elements in bec.psi are approximately 1
        condition = np.allclose(abs(bec.psi), 1, atol=tolerance)
        self.assertTrue(condition, "Elements in bec.psi are not approximately 1")


    # something is wrong with this. To fix
    # def test_vortex_tracker(self):
    #     bec = cf.BoseEinsteinCondensate(2,xRes=31,yRes=31)
    #     bec.conf_insert_vortex_dipole()
    #     bec.evolve_relax(200)
    #     bec.plot_field(abs(bec.psi))
    #     plt.show()
    #     Dnodes = bec.calc_vortex_nodes()
    #     print(Dnodes)
    #     self.assertEqual(len(Dnodes), 2)


if __name__ == '__main__':
    unittest.main()
