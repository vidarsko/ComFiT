import unittest
import sys
import os

# Adjust the path to import the comfit package
sys.path.append(os.path.abspath('../'))
import comfit as cf
import matplotlib.pyplot as plt
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

    ## 2D triangular tests

    def test_phase_field_crystal_2d_triangular_evolve_PFC_with_dislocation_dipole(self):
        """Test PhaseFieldCrystal2DTriangular calc_amplitudes_with_dislocation_dipole."""
        pfc = cf.PhaseFieldCrystal2DTriangular(21, 14)
        eta = pfc.calc_amplitudes_with_dislocation_dipole()
        pfc.conf_PFC_from_amplitudes(eta)
        pfc.evolve_PFC(200)

        # Check if pfc.psi is real-valued
        self.assertTrue(np.isrealobj(pfc.psi), "pfc.psi contains complex values.")

        # Check if pfc.psi does not contain NaN values
        self.assertFalse(np.any(np.isnan(pfc.psi)), "pfc.psi contains NaN values.")

    def test_phase_field_crystal_2d_triangular_initial_amplitudes(self):

        params = [
            {'nx': 1, 'ny': 1, 'r': -0.1, 't': 0, 'v': 1, 'psi0': -0.1},
            {'nx': 1, 'ny': 1, 'r': -0.2, 't': 0, 'v': 1, 'psi0': -0.2},
            {'nx': 1, 'ny': 1, 'r': -0.3, 't': 0, 'v': 1, 'psi0': -0.3},
            # Add more parameter combinations as needed
        ]

        # The following values were calculated using the legacy PFC matlab code from Vidars PhD
        As = [0.091180521680209,
              0.123266639978645,
              0.134833147735479]

        for p, A in zip(params, As):
            with self.subTest(p=p):
                try:
                    pfc = cf.PhaseFieldCrystal2DTriangular(**p)
                    pfc.conf_PFC_from_amplitudes()
                    self.assertAlmostEqual(pfc.A, A, places=5)
                except Exception as e:
                    self.fail(f"PFC Triangular amplitude calculation failed with {p}: {e}")

    def test_phase_field_crystal_2d_triangular_identify_dislocations(self):
        
        # Positions of the dislocations

        for dislocation_type in [1,2,3]:
            try:
                # Initialize PFC
                pfc = cf.PhaseFieldCrystal2DTriangular(21, 14)
                x1=pfc.xmax/3 
                y1=pfc.ymax/2
                x2=2*pfc.xmax/3
                y2=pfc.ymax/2

                eta = pfc.calc_amplitudes_with_dislocation_dipole(
                    dislocation_type=1,
                    x1=x1, y1=y1,
                    x2=x2, y2=y2)
                pfc.conf_PFC_from_amplitudes(eta)
                pfc.evolve_PFC(100)
                dislocation_nodes = pfc.calc_dislocation_nodes()

                # Check if there are two dislocations
                self.assertEqual(len(dislocation_nodes),2)
                
                # Check if the dislocations are correct
                for dislocation in dislocation_nodes:
                    self.assertEqual(np.linalg.norm(dislocation['Burgers_vector']), pfc.a0)
                    if (dislocation['Burgers_vector'] == -pfc.a[dislocation_type-1]).all():
                        self.assertAlmostEqual(dislocation['position'][0],x1,delta=pfc.a0)
                        self.assertAlmostEqual(dislocation['position'][1],y1,delta=pfc.a0)

                    elif (dislocation['Burgers_vector'] == pfc.a[dislocation_type-1]).all():
                        self.assertAlmostEqual(dislocation['position'][0],x2,delta=pfc.a0)
                        self.assertAlmostEqual(dislocation['position'][1],y2,delta=pfc.a0)

            except Exception as e:
                self.fail(f"PFC Triangular dislocation identification failed with dislocation type {dislocation_type}: {e}")

    ## 2D square tests

    def test_phase_field_crystal_2d_square_initial_amplitudes(self):

        params = [
            {'nx': 1, 'ny': 1, 'r': -0.1, 't': 0, 'v': 1, 'psi0': -0.1},
            {'nx': 1, 'ny': 1, 'r': -0.2, 't': 0, 'v': 1, 'psi0': -0.2},
            {'nx': 1, 'ny': 1, 'r': -0.3, 't': 0, 'v': 1, 'psi0': -0.3},
            # Add more parameter combinations as needed
        ]

        #The following values were calculated using the legacy PFC matlab code from Vidars PhD
        As = [0.098138558341561,
              0.132946527647919,
              0.152161548344077]

        Bs = [0.046917143364986,
              0.073844303366386,
              0.090505473985471]

        for p, A, B in zip(params, As, Bs):
            with self.subTest(p=p):
                try:
                    pfc = cf.PhaseFieldCrystal2DSquare(**p)
                    pfc.conf_PFC_from_amplitudes()
                    self.assertAlmostEqual(pfc.A, A, places=5)
                    self.assertAlmostEqual(pfc.B, B, places=5)
                except Exception as e:
                    self.fail(f"PFC Square amplitude calculation failed with {p}: {e}")

    def test_phase_field_crystal_2d_square_identify_dislocations(self):
        
        # Positions of the dislocations

        for dislocation_type in [1,2]:
            try:
                # Initialize PFC
                pfc = cf.PhaseFieldCrystal2DSquare(21, 21)
                x1=pfc.xmax/3 
                y1=pfc.ymax/2
                x2=2*pfc.xmax/3
                y2=pfc.ymax/2

                eta = pfc.calc_amplitudes_with_dislocation_dipole(
                    dislocation_type=1,
                    x1=x1, y1=y1,
                    x2=x2, y2=y2)
                pfc.conf_PFC_from_amplitudes(eta)
                pfc.evolve_PFC(100)
                dislocation_nodes = pfc.calc_dislocation_nodes()

                # Check if there are two dislocations
                self.assertEqual(len(dislocation_nodes),2)
                
                # Check if the dislocations are correct
                for dislocation in dislocation_nodes:
                    self.assertEqual(np.linalg.norm(dislocation['Burgers_vector']), pfc.a0)
                    if (dislocation['Burgers_vector'] == -pfc.a[dislocation_type-1]).all():
                        self.assertAlmostEqual(dislocation['position'][0],x1,delta=pfc.a0)
                        self.assertAlmostEqual(dislocation['position'][1],y1,delta=pfc.a0)

                    elif (dislocation['Burgers_vector'] == pfc.a[dislocation_type-1]).all():
                        self.assertAlmostEqual(dislocation['position'][0],x2,delta=pfc.a0)
                        self.assertAlmostEqual(dislocation['position'][1],y2,delta=pfc.a0)

            except Exception as e:
                self.fail(f"PFC Square dislocation identification failed with dislocation type {dislocation_type}: {e}")


    def test_phase_field_crystal_3d_body_centered_cubic_initial_amplitudes(self):

        params = [
            {'nx': 1, 'ny': 1, 'nz': 1, 'r': -0.1, 't': 0, 'v': 1, 'psi0': -0.125},
            {'nx': 1, 'ny': 1, 'nz': 1, 'r': -0.2, 't': 0, 'v': 1, 'psi0': -0.225},
            {'nx': 1, 'ny': 1, 'nz': 1, 'r': -0.3, 't': 0, 'v': 1, 'psi0': -0.325},
            # Add more parameter combinations as needed
        ]

        # The following values were calculated using the legacy PFC matlab code from Vidars PhD
        As = [0.054854797457965,
              0.074378423185648,
              0.082099011165377]

        for p, A in zip(params, As):
            with self.subTest(p=p):
                try:
                    pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(**p)
                    pfc.conf_PFC_from_amplitudes()
                    self.assertAlmostEqual(pfc.A, A, places=5)
                except Exception as e:
                    self.fail(f"PFC BCC amplitude calculation failed with {p}: {e}")

    


    def test_phase_field_crystal_3d_face_centered_cubic_initial_amplitudes(self):

        params = [
            {'nx': 1, 'ny': 1, 'nz': 1, 'r': -0.1, 't': 0, 'v': 1, 'psi0': -0.125},
            {'nx': 1, 'ny': 1, 'nz': 1, 'r': -0.2, 't': 0, 'v': 1, 'psi0': -0.225},
            {'nx': 1, 'ny': 1, 'nz': 1, 'r': -0.3, 't': 0, 'v': 1, 'psi0': -0.325},
            # Add more parameter combinations as needed
        ]

        # The following values were calculated using the legacy PFC matlab code from Vidars PhD
        As = [0.056481938628754,
              0.077467125827116,
              0.088828545091084]

        Bs = [0.038972621699196,
              0.056347683829042,
              0.066589290566421]

        for p, A, B in zip(params, As, Bs):
            with self.subTest(p=p):
                try:
                    pfc = cf.PhaseFieldCrystal3DFaceCenteredCubic(**p)
                    pfc.conf_PFC_from_amplitudes()
                    self.assertAlmostEqual(pfc.A, A, places=5)
                    self.assertAlmostEqual(pfc.B, B, places=5)
                except Exception as e:
                    self.fail(f"PFC FCC amplitude calculation failed with {p}: {e}")

    def test_phase_field_crystal_3d_body_centered_cubic_identify_dislocations(self):
        # Set the plane equation
        a = 0.3
        b = 0.7
        c = 0.1
        
        for dislocation_type in range(1,5):
            try:
                pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(nx=13, ny=13, nz=13)
                # Insert a dislocation ring
                eta = pfc.calc_amplitudes_with_dislocation_ring(
                    normal_vector=[a,b,c],
                    position=pfc.rmid,
                    radius=pfc.xmax/3,
                    dislocation_type=dislocation_type
                )
                pfc.conf_PFC_from_amplitudes(eta)

                # Relax the system
                pfc.evolve_PFC(20)
                dislocation_nodes = pfc.calc_dislocation_nodes()
                
                #Check that there are nodes
                self.assertGreater(len(dislocation_nodes),0,"Dislocation nodes not found")

                #Check that all the nodes belong to the correct plane
                for dislocation_node in dislocation_nodes:
                    # Calculate the plane equation
                    plane_equation = a*(dislocation_node['position'][0]-pfc.rmid[0]) \
                                    + b*(dislocation_node['position'][1]-pfc.rmid[1]) \
                                    + c*(dislocation_node['position'][2]-pfc.rmid[2])

                    # Check if the plane equation is zero
                    self.assertAlmostEqual(plane_equation,0,delta=0.5*pfc.a0)
                    
                    # Check that the Burgers vector is correct
                    print(dislocation_node['Burgers_vector'])
                    print(pfc.a[dislocation_type-1])
                    self.assertTrue((dislocation_node['Burgers_vector'] == pfc.a[dislocation_type-1]).all())

            except Exception as e:
                self.fail(f"PFC BCC dislocation identification failed with dislocation type {dislocation_type}: {e}")

    def test_phase_field_crystal_3d_face_centered_cubic_initial_amplitudes(self):

        params = [
            {'nx': 1, 'ny': 1, 'nz': 1, 'r': -0.1, 't': 0, 'v': 1, 'psi0': -0.125},
            {'nx': 1, 'ny': 1, 'nz': 1, 'r': -0.2, 't': 0, 'v': 1, 'psi0': -0.225},
            {'nx': 1, 'ny': 1, 'nz': 1, 'r': -0.3, 't': 0, 'v': 1, 'psi0': -0.325},
            # Add more parameter combinations as needed
        ]

        # The following values were calculated using the legacy PFC matlab code from Vidars PhD
        As = [0.072799451159192,
              0.082792401061476,
              0.090170308517259]

        Bs = [0.020831100850532,
              0.037148906733000,
              0.047858014941568]

        Cs = [-0.004840409953021,
              0.011422779247214,
              0.022268287314380]

        for p, A, B, C in zip(params, As, Bs, Cs):
            with self.subTest(p=p):
                try:
                    pfc = cf.PhaseFieldCrystal3DSimpleCubic(**p)
                    pfc.conf_PFC_from_amplitudes()
                    self.assertAlmostEqual(pfc.A, A, places=5)
                    self.assertAlmostEqual(pfc.B, B, places=5)
                    self.assertAlmostEqual(pfc.C, C, places=5)
                except Exception as e:
                    self.fail(f"PFC Simple Cubic amplitude calculation failed with {p}: {e}")


if __name__ == '__main__':
    unittest.main()
    
