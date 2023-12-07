import unittest
import sys
import os

# Adjust the path to import the comfit package
sys.path.append(os.path.abspath('../'))
import comfit as cf

class TestBoseEinsteinCondensate(unittest.TestCase):

    def test_init_with_dimension(self):
        """Test BoseEinsteinCondensate initialization with a dimension parameter."""
        for dim in [1, 2, 3]:
            try:
                bec = cf.BoseEinsteinCondensate(dim)
                self.assertIsInstance(bec, cf.BoseEinsteinCondensate)
            except Exception as e:
                self.fail(f"Initialization failed with dimension {dim}: {e}")

if __name__ == '__main__':
    unittest.main()
