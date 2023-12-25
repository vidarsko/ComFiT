import unittest
import sys
import os

# Adjust the path to import the comfit package
sys.path.append(os.path.abspath('../'))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

class TestNematicLiquidCrystal(unittest.TestCase):

    def test_init_with_dimension(self):
        """Test NematicLiquidCrystal initialization with a dimension parameter."""
        for dim in [1, 2, 3]:
            try:
                nlc = cf.NematicLiquidCrystal(dim)
                self.assertIsInstance(nlc, cf.NematicLiquidCrystal)
            except Exception as e:
                self.fail(f"Initialization failed with dimension {dim}: {e}")

    def test_no_flow_evolver(self):
        """Test the enm.evolve_nematic_no_flow"""
        nem = cf.NematicLiquidCrystal(2, xRes=13, yRes=4)
        nem.conf_initial_condition_disordered()
        nem.evolve_nematic_no_flow(300)

        # Set the tolerance for approximation
        tolerance = 0.01

        # Check if all elements in bec.psi are approximately 1
        S = np.sqrt(nem.B)/2
        condition = np.allclose(abs(nem.Q[0][0]),S , atol=tolerance)
        self.assertTrue(condition, "Elements in nem.Q are not approximately 1")


if __name__ == '__main__':
    unittest.main()
