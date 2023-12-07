import unittest
import sys

# Adjust the path to import your comfit package
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

if __name__ == '__main__':
    unittest.main()
