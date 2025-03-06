import unittest
import sys
import os

# Run 
# Adjust the path to import the comfit package
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
import comfit as cf
import numpy as np
import scipy as sp

class TestBaseSystem(unittest.TestCase):

    def tearDown(self):
        import matplotlib.pyplot as plt
        plt.close('all')

    def test_plot_field_no_parameters(self):
        """ Test plotting of field """

        # 1 dimension
        # Initialize BaseSystem object
        bs = cf.BaseSystem(1, xRes=32, dx=0.1)

        # Create field
        np.random.seed(20423536)
        field = np.random.rand(bs.xRes)

        # Plot field
        try:
            bs.plot_field(field)
            
        except Exception as e:
            self.fail(f"Plotting failed: {e}")

        # 2 dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(2, xRes=32, dx=0.1, yRes=32, dy=0.1)

        # Create field
        np.random.seed(76427622)
        field = np.random.rand(bs.yRes, bs.xRes)

        # Plot field
        try:
            bs.plot_field(field)
            
        except Exception as e:
            self.fail(f"Plotting failed: {e}")

        # 3 dimensions
        # Initialize BaseSystem object

        bs = cf.BaseSystem(3, xRes=32, dx=0.1, yRes=32, dy=0.1, zRes=32, dz=0.1)

        # Create field
        np.random.seed(42000637)
        field = np.random.rand(bs.zRes, bs.yRes, bs.xRes)

        # Plot field
        try:
            bs.plot_field(field)
            
        except Exception as e:
            self.fail(f"Plotting failed: {e}")

    def test_plot_complex_field_no_parameters(self):
        """ Test plotting of complex field """

        # 1 dimension
        # Initialize BaseSystem object
        bs = cf.BaseSystem(1, xRes=32, dx=0.1)

        # Create field
        np.random.seed(83618706)
        field = np.random.rand(bs.xRes) + 1j*np.random.rand(bs.xRes)

        # Plot field
        try:
            bs.plot_complex_field(field)
            
        except Exception as e:
            self.fail(f"Plotting failed: {e}")

        # 2 dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(2, xRes=32, dx=0.1, yRes=32, dy=0.1)

        # Create field
        np.random.seed(86050252)
        field = np.random.rand(bs.yRes, bs.xRes) + 1j*np.random.rand(bs.yRes, bs.xRes)

        # Plot field
        try:
            bs.plot_complex_field(field)
            
        except Exception as e:
            self.fail(f"Plotting failed: {e}")

        # 3 dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(3, xRes=32, dx=0.1, yRes=32, dy=0.1, zRes=32, dz=0.1)

        # Create field
        np.random.seed(63363579)
        field = np.random.rand(bs.zRes, bs.yRes, bs.xRes) + 1j*np.random.rand(bs.zRes, bs.yRes, bs.xRes)

        # Plot field
        try:
            bs.plot_complex_field(field)
            
        except Exception as e:
            self.fail(f"Plotting failed: {e}")
    
    def test_plot_angle_field_no_parameters(self):
        """ Test plotting of angle field """

        # 1 dimension
        # Initialize BaseSystem object
        bs = cf.BaseSystem(1, xRes=32, dx=0.1)

        # Calculate angle field
        angle_field = 2*np.pi*np.random.rand(bs.xRes)-np.pi

        # Plot field
        try:
            bs.plot_angle_field(angle_field)
            
        except Exception as e:
            self.fail(f"Plotting failed: {e}")


        # 2 dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(2, xRes=32, dx=0.1, yRes=32, dy=0.1)

        # Calculate angle field
        np.random.seed(58271391)
        angle_field = 2*np.pi*np.random.rand(bs.yRes, bs.xRes)-np.pi

        # Plot field
        try:
            bs.plot_angle_field(angle_field)
            
        except Exception as e:
            self.fail(f"Plotting failed: {e}")

        # 3 dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(3, xRes=32, dx=0.1, yRes=32, dy=0.1, zRes=32, dz=0.1)

        # Calculate angle field
        np.random.seed(77579665)
        angle_field = 2*np.pi*np.random.rand(bs.zRes, bs.yRes, bs.xRes)-np.pi

        # Plot field
        try:
            bs.plot_angle_field(angle_field)
            
        except Exception as e:
            self.fail(f"Plotting failed: {e}")

    def test_plot_vector_field_no_parameters(self):
        """ Test plotting of vector field """

        # 1 dimension
        # Initialize BaseSystem object
        bs = cf.BaseSystem(1, xRes=32, dx=0.1)

        for vector_dim in range(1,4):
        # Create field
            np.random.seed(52137189)
            field = np.random.rand(vector_dim, bs.xRes)

            # Plot field
            try:
                bs.plot_vector_field(field)
                
            except Exception as e:
                self.fail(f"Plotting failed: {e}")

        # 2 dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(2, xRes=32, dx=0.1, yRes=32, dy=0.1)

        for vector_dim in range(1,4):
            # Create field
            np.random.seed(12138881)
            field = np.random.rand(vector_dim, bs.yRes, bs.xRes)

            # Plot field
            try:
                bs.plot_vector_field(field)
                
            except Exception as e:
                self.fail(f"Plotting failed: {e}")

        # 3 dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(3, xRes=32, dx=0.1, yRes=32, dy=0.1, zRes=32, dz=0.1)


        for vector_dim in range(1,4):
            # Create field
            np.random.seed(19623034)
            field = np.random.rand(vector_dim, bs.zRes, bs.yRes, bs.xRes)

            # Plot field
            try:
                bs.plot_vector_field(field)
                
            except Exception as e:
                self.fail(f"Plotting failed: {e}")

    def test_plot_field_in_plane(self):

        # 3 dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(3, xRes=32, dx=0.1, yRes=32, dy=0.1, zRes=32, dz=0.1)

        # Create field
        np.random.seed(14803969)
        field = np.random.rand(bs.zRes, bs.yRes, bs.xRes)

        # Plot field
        try:
            bs.plot_field_in_plane(field)
        except Exception as e:
            self.fail(f"Plotting failed: {e}")

    def test_plot_complex_field_in_plane(self):
        """ Test plotting of complex field in plane """

        # 3 dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(3, xRes=32, dx=0.1, yRes=32, dy=0.1, zRes=32, dz=0.1)

        # Create field
        np.random.seed(73158863)
        field = np.random.rand(bs.zRes, bs.yRes, bs.xRes) + 1j*np.random.rand(bs.zRes, bs.yRes, bs.xRes)

        # Plot field
        try:
            bs.plot_complex_field_in_plane(field)
        except Exception as e:
            self.fail(f"Plotting failed: {e}")

    def plot_angle_field_in_plane(self):
        """ Test plotting of angle field in plane """

        # 3 dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(3, xRes=32, dx=0.1, yRes=32, dy=0.1, zRes=32, dz=0.1)

        # Calculate angle field
        np.random.seed(15735094)
        angle_field = 2*np.pi*np.random.rand(bs.zRes, bs.yRes, bs.xRes)-np.pi

        # Plot field
        try:
            bs.plot_angle_field_in_plane(angle_field)
            
        except Exception as e:
            self.fail(f"Plotting failed: {e}")

    def plot_vector_field_in_plane(self):
        """ Test plotting of vector field in plane """

        # 3 dimensions
        # Initialize BaseSystem object
        bs = cf.BaseSystem(3, xRes=32, dx=0.1, yRes=32, dy=0.1, zRes=32, dz=0.1)

        for vector_dim in range(1,4):
            # Create field
            np.random.seed(56688066)
            field = np.random.rand(vector_dim, bs.zRes, bs.yRes, bs.xRes)

            # Plot field
            try:
                bs.plot_vector_field_in_plane(field)
                
            except Exception as e:
                self.fail(f"Plotting failed: {e}")

    def test_plot_field_with_different_parameters(self):
        bs = cf.BaseSystem(3, xRes=32, dx=0.1, yRes=32, dy=0.1, zRes=32, dz=0.1)
        bs.a0 = 2

        # Create field
        np.random.seed(97623584)
        field = np.random.rand(bs.zRes, bs.yRes, bs.xRes)

        # Test axis labels
        fig, ax = bs.plot_field(field, xlabel='x', ylabel='y', zlabel='z')
        self.assertEqual(fig['layout']['scene']['xaxis']['title']['text'], 'x')
        self.assertEqual(fig['layout']['scene']['yaxis']['title']['text'], 'y')
        self.assertEqual(fig['layout']['scene']['zaxis']['title']['text'], 'z')

        # Test xmin, xmax, ymin, ymax, zmin, zmax
        fig, ax = bs.plot_field(field, xmin=-3, xmax=4, ymin=-5, ymax=6, zmin=-7, zmax=8)
        self.assertEqual(fig['layout']['scene']['xaxis']['range'], (-3/bs.a0, 4/bs.a0))
        self.assertEqual(fig['layout']['scene']['yaxis']['range'], (-5/bs.a0, 6/bs.a0))
        self.assertEqual(fig['layout']['scene']['zaxis']['range'], (-7/bs.a0, 8/bs.a0))

        # Test xlim, ylim, zlim
        fig, ax = bs.plot_field(field, xlim=(-3, 4), ylim=(-5, 6), zlim=(-7, 8))
        self.assertEqual(fig['layout']['scene']['xaxis']['range'], (-3/bs.a0, 4/bs.a0))
        self.assertEqual(fig['layout']['scene']['yaxis']['range'], (-5/bs.a0, 6/bs.a0))
        self.assertEqual(fig['layout']['scene']['zaxis']['range'], (-7/bs.a0, 8/bs.a0))

        # Test xticks, yticks, zticks
        fig, ax = bs.plot_field(field, xticks=[-3, 0, 4], yticks=[-5, 0, 6], zticks=[-7, 0, 8])
        self.assertEqual(fig['layout']['scene']['xaxis']['tickvals'], (-3, 0, 4))
        self.assertEqual(fig['layout']['scene']['yaxis']['tickvals'], (-5, 0, 6))
        self.assertEqual(fig['layout']['scene']['zaxis']['tickvals'], (-7, 0, 8))


if __name__ == '__main__':
    unittest.main()
