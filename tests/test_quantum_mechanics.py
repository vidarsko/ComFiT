import unittest
import sys
import os

# Adjust the path to import the comfit package
sys.path.append(os.path.abspath('../'))
import comfit as cf


class TestQuantumMechanics(unittest.TestCase):

    def test_init_with_dimension(self):
        """Test QuantumMechanics initialization with a dimension parameter."""
        for dim in [1, 2, 3]:
            try:
                qm = cf.QuantumMechanics(dim)
                self.assertIsInstance(qm, cf.QuantumMechanics)
            except Exception as e:
                self.fail(f"Initialization failed with dimension {dim}: {e}")

        # Make a test to check that the QM mass is conserved

if __name__ == '__main__':
    unittest.main()
