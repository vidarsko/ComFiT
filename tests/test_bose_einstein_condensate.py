import unittest
import sys
import os

# Adjust the path to import the comfit package
sys.path.append(os.path.abspath('../'))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

class TestBoseEinsteinCondensate(unittest.TestCase):
    def test_init_with_dimension(self):
        """Test BoseEinsteinCondensate initialization with a dimension parameter."""
        for dim in [1, 2, 3]:
            try:
                bec = cf.BoseEinsteinCondensate(dim)
                self.assertIsInstance(bec, cf.BoseEinsteinCondensate)
            except Exception as e:
                self.fail(f"Initialization failed with dimension {dim}: {e}")

    def test_dGPE_relaxer(self):
        """Test the dGPE relaxer."""
        bec = cf.BoseEinsteinCondensate(2,xRes=13,yRes=4)
        bec.conf_initial_condition_disordered()
        bec.evolve_relax(300)

        # Set the tolerance for approximation
        tolerance = 0.01

        # Check if all elements in bec.psi are approximately 1
        condition = np.allclose(abs(bec.psi), 1, atol=tolerance)
        self.assertTrue(condition, "Elements in bec.psi are not approximately 1")

    def test_vortex_tracker(self):
        bec = cf.BoseEinsteinCondensate(2,xRes=31,yRes=31)
        bec.conf_insert_vortex_dipole(
            dipole_vector=[bec.xmax/3,0],
            dipole_position=bec.rmid)
        bec.evolve_relax(100)
        Dnodes = bec.calc_vortex_nodes()

        self.assertEqual(len(Dnodes),2,"Dipole nodes not found")
        xpos = [node['position'][0] for node in Dnodes]
        ypos = [node['position'][1] for node in Dnodes]

        # Sort the lists according to xpos
        sorted_lists = sorted(zip(xpos, ypos))

        # Unzip the sorted lists
        xpos, ypos = zip(*sorted_lists)
        
        # Testing the position of the first node
        self.assertAlmostEqual(xpos[0],bec.xmax/3,delta=1)
        self.assertAlmostEqual(ypos[0],bec.ymax/2,delta=1)

        # Testing the position of the second node
        self.assertAlmostEqual(xpos[1],2*bec.xmax/3,delta=1)
        self.assertAlmostEqual(ypos[1],bec.ymax/2,delta=1)







if __name__ == '__main__':
    unittest.main()
