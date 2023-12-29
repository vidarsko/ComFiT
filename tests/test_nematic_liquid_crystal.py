import unittest
import sys
import os

# Adjust the path to import the comfit package
sys.path.append(os.path.abspath('../'))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

class TestNematicLiquidCrystal(unittest.TestCase):

    def test_init_with_dimension(self):
        """Test NematicLiquidCrystal initialization with a dimension parameter."""
        for dim in [1, 2, 3]:
            try:
                nlc = cf.NematicLiquidCrystal(dim)
                self.assertIsInstance(nlc, cf.NematicLiquidCrystal)
            except Exception as e:
                self.fail(f"Initialization failed with dimension {dim}: {e}")

    def test_no_flow_evolver(self):
        """Test the enm.evolve_nematic_no_flow"""
        nem = cf.NematicLiquidCrystal(2, xRes=13, yRes=4)
        nem.conf_initial_condition_disordered()
        nem.evolve_nematic_no_flow(300)

        # Set the tolerance for approximation
        tolerance = 0.01

        # Check if all elements in bec.psi are approximately 1
        S = np.sqrt(nem.B)/2
        condition = np.allclose(abs(nem.Q[0]),S , atol=tolerance)
        self.assertTrue(condition, "Elements in nem.Q are not approximately 1")


    def test_nematic_evolver_with_defect_dipole(self):
        """Test the enm.evolve_nematic_no_flow"""
        nem = cf.NematicLiquidCrystal(2, xRes=31, yRes=31)
        nem.conf_insert_vortex_dipole(dipole_vector=[nem.xmax/3,0],
            dipole_position=nem.rmid)
        nem.evolve_nematic_no_flow(10)
        nem.evolve_nematic(1)
        Dnodes = nem.calc_vortex_nodes_nem()

        self.assertEqual(len(Dnodes), 2, "Dipole nodes not found")

        # Create lists of the x and y positions of the nodes
        xpos = [node['position'][0] for node in Dnodes]
        ypos = [node['position'][1] for node in Dnodes]

        # Sort the lists according to xpos
        sorted_lists = sorted(zip(xpos, ypos))

        # Unzip the sorted lists
        xpos, ypos = zip(*sorted_lists)

        # Testing the position of the first node
        self.assertAlmostEqual(xpos[0], nem.xmax / 3, delta=1)
        self.assertAlmostEqual(ypos[0], nem.ymax / 2, delta=1)

        # Testing the position of the second node
        self.assertAlmostEqual(xpos[1], 2 * nem.xmax / 3, delta=1)
        self.assertAlmostEqual(ypos[1], nem.ymax / 2, delta=1)




if __name__ == '__main__':
    unittest.main()
