import unittest
import sys
import os

# Run 
# Adjust the path to import the comfit package
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np


class TestQuantumMechanics(unittest.TestCase):

    def test_init_with_dimension(self):
        """Test QuantumMechanics initialization with a dimension parameter."""
        for dim in [1, 2, 3]:
            try:
                qm = cf.QuantumMechanics(dim)
                self.assertIsInstance(qm, cf.QuantumMechanics)
            except Exception as e:
                self.fail(f"Initialization failed with dimension {dim}: {e}")

    
    def test_evolution_conserved(self):
        """ Test that the evolution of a quantum state conserves the norm."""

        params = [{},{},{'xRes': 30, 'yRes': 30, 'zRes': 30}]

        for dim, p in zip([1, 2, 3],params):
            
            # Initialize a quantum mechanics system
            qm = cf.QuantumMechanics(dim,**p)
            qm.conf_initial_condition_Gaussian()

            # qm.plot_complex_field(qm.psi)
            # plt.show()

            # Check that the norm is approximately 1 at beginning
            norm_at_time0 = qm.calc_integrate_field(abs(qm.psi)**2)
            self.assertAlmostEqual(norm_at_time0, 1.0, places=3)
            
            # Evolve the state and check that the norm is still approximately 1
            qm.evolve_schrodinger(1000)
            norm_at_time1 = qm.calc_integrate_field(abs(qm.psi)**2)
            self.assertAlmostEqual(norm_at_time1, 1.0, places=3)
        


if __name__ == '__main__':
    unittest.main()
