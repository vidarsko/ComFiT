import unittest
import sys
import os

# Adjust the path to import the comfit package
sys.path.append(os.path.abspath('../'))
import comfit as cf

import numpy as np

class TestPhaseFieldCrystal(unittest.TestCase):

    def test_phase_field_crystal_2d_triangular_init(self):
        """Test PhaseFieldCrystal2DTriangular initialization with nx=ny=1."""
        try:
            pfc1d = cf.PhaseFieldCrystal2DTriangular(1, 1)
            self.assertIsInstance(pfc1d, cf.PhaseFieldCrystal2DTriangular)
        except Exception as e:
            self.fail(f"Initialization of PhaseFieldCrystal1DPeriodic failed: {e}")

    def test_phase_field_crystal_2d_triangular_init_with_params(self):
        """Test PhaseFieldCrystal2DTriangular initialization with different parameters."""
        params = [
            {'nx': 11, 'ny': 14, 'dt': 0.1},
            {'nx': 32, 'ny': 32, 'dt': 0.3},
            # Add more parameter combinations as needed
        ]
        for p in params:
            with self.subTest(p=p):
                try:
                    pfc = cf.PhaseFieldCrystal2DTriangular(**p)
                    self.assertIsInstance(pfc, cf.PhaseFieldCrystal2DTriangular)
                except Exception as e:
                    self.fail(f"Initialization failed with parameters {p}: {e}")

    # Add similar test methods for other subclasses like PhaseFieldCrystal2DTriangular, etc.

    def test_phase_field_crystal_2d_triangular_evolve_PFC_with_dislocation_dipole(self):
        """Test PhaseFieldCrystal2DTriangular calc_amplitudes_with_dislocation_dipole."""
        pfc = cf.PhaseFieldCrystal2DTriangular(21,14)
        eta = pfc.calc_amplitudes_with_dislocation_dipole()
        pfc.conf_PFC_from_amplitudes(eta)
        pfc.evolve_PFC(200)

        # Check if pfc.psi is real-valued
        self.assertTrue(np.isrealobj(pfc.psi), "pfc.psi contains complex values.")

        # Check if pfc.psi does not contain NaN values
        self.assertFalse(np.any(np.isnan(pfc.psi)), "pfc.psi contains NaN values.")

if __name__ == '__main__':
    unittest.main()
