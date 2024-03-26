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
        np.random.seed(29618953)
        nem.conf_initial_condition_ordered(noise_strength=2)
        nem.evolve_nematic_no_flow(500)

        # Set the tolerance for approximation
        tolerance = 0.01

        # Check if nem.Q have approacjed the expected=
        S_0 = np.sqrt(nem.B)
        S,n = nem.calc_order_and_director()
        condition = np.allclose(S_0,S , atol=tolerance)
        self.assertTrue(condition, "Elements in nem.Q are not relaxed to the correct value")

    def test_no_flow_evolver_3D(self):
        for i in range(10):
            c = 2*np.random.rand()
            nem = cf.NematicLiquidCrystal(3, xRes=13, yRes=4,zRes=13,C = c)
            np.random.seed(29820894)
            nem.conf_initial_condition_ordered(noise_strength=2)

            nem.evolve_nematic_no_flow(300)


            # Set the tolerance for approximation
            tolerance = 0.01

            # Check if all elements in  are approximately 1
            S_0 =1/8* nem.C/nem.A + 1/2 * np.sqrt(nem.C**2 /(16*nem.A**2) + 3*nem.B)

            S, n = nem.calc_order_and_director()
            condition = np.allclose(S_0, S, atol=tolerance)
            self.assertTrue(condition, "Elements in nem.Q are not relaxed to the correct value")



    def test_nematic_evolver_with_defect_dipole(self):
        """Test the enm.evolve_nematic_no_flow"""
        nem = cf.NematicLiquidCrystal(2, xRes=31, yRes=31)
        nem.conf_insert_disclination_dipole(dipole_vector=[nem.xmax/3,0],
            dipole_position=nem.rmid)
        nem.evolve_nematic_no_flow(10)
        nem.evolve_nematic(1)
        Dnodes = nem.calc_disclination_nodes_nem()

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
