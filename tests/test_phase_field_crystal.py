import unittest
import sys
import os

# Adjust the path to import the comfit package
sys.path.append(os.path.abspath('../'))
import comfit as cf

class TestPhaseFieldCrystal(unittest.TestCase):

    def test_phase_field_crystal_init(self):
        """Test PhaseFieldCrystal initialization with a dimension parameter."""
        for dim in [1, 2, 3]:
            try:
                pfc = cf.PhaseFieldCrystal(dimension=dim)
                self.assertIsInstance(pfc, cf.PhaseFieldCrystal)
            except Exception as e:
                self.fail(f"Initialization of PhaseFieldCrystal failed with dimension {dim}: {e}")

    def test_phase_field_crystal_2d_triangular_init(self):
        """Test PhaseFieldCrystal2DTriangular initialization with nx=ny=1."""
        try:
            pfc1d = cf.PhaseFieldCrystal1DPeriodic(1, 1)
            self.assertIsInstance(pfc1d, cf.PhaseFieldCrystal2DTriangular)
        except Exception as e:
            self.fail(f"Initialization of PhaseFieldCrystal1DPeriodic failed: {e}")

    # Add similar test methods for other subclasses like PhaseFieldCrystal2DTriangular, etc.

if __name__ == '__main__':
    unittest.main()
