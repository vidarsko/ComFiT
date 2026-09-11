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
            qm.conf_initial_condition_gaussian()

            # qm.plot_complex_field(qm.psi)
            # plt.show()

            # Check that the norm is approximately 1 at beginning
            norm_at_time0 = qm.calc_integrate_field(abs(qm.psi)**2)
            self.assertAlmostEqual(norm_at_time0, 1.0, places=3)
            
            # Evolve the state and check that the norm is still approximately 1
            qm.evolve_schrodinger(1000)
            norm_at_time1 = qm.calc_integrate_field(abs(qm.psi)**2)
            self.assertAlmostEqual(norm_at_time1, 1.0, places=3)

    def test_psi_component_axis(self):
        """qm.psi always carries a leading component axis of size 1."""
        qm = cf.QuantumMechanics(2, xRes=21, yRes=21)

        qm.conf_initial_condition_gaussian()
        self.assertEqual(qm.psi.shape, (1, qm.xRes, qm.yRes))
        self.assertEqual(qm.psi_f.shape, qm.psi.shape)

        qm.conf_wavefunction(qm.psi[0])
        self.assertEqual(qm.psi.shape, (1, qm.xRes, qm.yRes))
        self.assertEqual(qm.psi_f.shape, qm.psi.shape)


if __name__ == '__main__':
    unittest.main()
