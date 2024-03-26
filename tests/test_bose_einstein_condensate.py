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
        bec = cf.BoseEinsteinCondensate(2,xRes=23,yRes=14)
        np.random.seed(29160566)
        bec.conf_initial_condition_disordered()
        bec.evolve_relax(1000)

        # Set the tolerance for approximation
        tolerance = 0.01

        # Check if all elements in bec.psi are approximately 1
        condition = np.allclose(abs(bec.psi), 1, atol=tolerance)
        self.assertTrue(condition, "Elements in bec.psi are not approximately 1")

    def test_vortex_tracker_2D(self):
        bec = cf.BoseEinsteinCondensate(2,xRes=31,yRes=31)
        bec.conf_insert_vortex_dipole(
            dipole_vector=[bec.xmax/3,0],
            dipole_position=bec.rmid)
        bec.evolve_relax(100)
        Dnodes = bec.calc_vortex_nodes()

        # Testing the number of nodes
        self.assertEqual(len(Dnodes),2,"Dipole nodes not found")

        # Create lists of the x and y positions of the nodes
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

    def test_vortex_tracker_3D(self):
        bec = cf.BoseEinsteinCondensate(3,xRes=31, yRes=31, zRes=31)
        # Set the plane equation
        a = 0.2
        b = 0.3
        c = 0.5
        
        # Insert a vortex ring
        bec.conf_insert_vortex_ring(
            normal_vector=[a,b,c],
            position=bec.rmid, 
            radius=bec.xmax/3)
        
        # Relax the system
        bec.evolve_relax(100)
        Dnodes = bec.calc_vortex_nodes()
        
        #Check that there are nodes
        self.assertGreater(len(Dnodes),0,"Vortex nodes not found")

        bec.plot_vortex_nodes(Dnodes)
        plt.show()
        #Check that all the nodes belong to the correct plane
        for node in Dnodes:
            # Calculate the plane equation
            plane_equation = a*(node['position'][0]-bec.rmid[0]) \
                            + b*(node['position'][1]-bec.rmid[1]) \
                            + c*(node['position'][2]-bec.rmid[2])

            # Check if the plane equation is zero
            self.assertAlmostEqual(plane_equation,0,delta=0.5)
            







if __name__ == '__main__':
    unittest.main()
