import unittest
import sys

# Adjust the path to import your comfit package
sys.path.append(os.path.abspath('../'))
import comfit as cf


class TestNematicLiquidCrystal(unittest.TestCase):

    def test_init_with_dimension(self):
        """Test NematicLiquidCrystal initialization with a dimension parameter."""
        for dim in [1, 2, 3]:
            try:
                nlc = cf.NematicLiquidCrystal(dim)
                self.assertIsInstance(nlc, cf.NematicLiquidCrystal)
            except Exception as e:
                self.fail(f"Initialization failed with dimension {dim}: {e}")

if __name__ == '__main__':
    unittest.main()
