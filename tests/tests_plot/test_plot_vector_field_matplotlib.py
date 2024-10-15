import unittest
import sys
import os
import numpy as np
import matplotlib

# Adjust the path to import the comfit package
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
import comfit as cf

class TestPlotVectorFieldMatplotlib(unittest.TestCase):

    def setUp(self):
        """Initialize the BaseSystem object for testing."""
        self.system1D = cf.BaseSystem(1, xRes=20)
        self.system2D = cf.BaseSystem(2, xRes=20, yRes=20)
        self.system3D = cf.BaseSystem(3, xRes=10, yRes=10, zRes=10)

        self.system1D.plot_lib = 'matplotlib'
        self.system2D.plot_lib = 'matplotlib'
        self.system3D.plot_lib = 'matplotlib'

    def test_plot_1D_vector_dim1(self):
        """Test 1D vector field with vector dimension 1."""
        vector_field = (np.random.rand(self.system1D.xRes),)
        fig, ax = self.system1D.plot_vector_field(vector_field)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_plot_1D_vector_dim2(self):
        """Test 1D vector field with vector dimension 2."""
        vector_field = (np.random.rand(self.system1D.xRes), np.random.rand(self.system1D.xRes))
        fig, ax = self.system1D.plot_vector_field(vector_field)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_plot_1D_vector_dim3(self):
        """Test 1D vector field with vector dimension 3."""
        vector_field = (
            np.random.rand(self.system1D.xRes),
            np.random.rand(self.system1D.xRes),
            np.random.rand(self.system1D.xRes)
        )
        fig, ax = self.system1D.plot_vector_field(vector_field)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_plot_2D_vector_dim1(self):
        """Test 2D vector field with vector dimension 1."""
        vector_field = (np.random.rand(self.system2D.xRes, self.system2D.yRes),)
        fig, ax = self.system2D.plot_vector_field(vector_field)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_plot_2D_vector_dim2(self):
        """Test 2D vector field with vector dimension 2."""
        vector_field = (
            np.random.rand(self.system2D.xRes, self.system2D.yRes),
            np.random.rand(self.system2D.xRes, self.system2D.yRes)
        )
        fig, ax = self.system2D.plot_vector_field(vector_field)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_plot_2D_vector_dim3(self):
        """Test 2D vector field with vector dimension 3."""
        vector_field = (
            np.random.rand(self.system2D.xRes, self.system2D.yRes),
            np.random.rand(self.system2D.xRes, self.system2D.yRes),
            np.random.rand(self.system2D.xRes, self.system2D.yRes)
        )
        fig, ax = self.system2D.plot_vector_field(vector_field)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_plot_3D_vector_dim1(self):
        """Test 3D vector field with vector dimension 1."""
        vector_field = (np.random.rand(self.system3D.xRes, self.system3D.yRes, self.system3D.zRes),)
        fig, ax = self.system3D.plot_vector_field(vector_field)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_plot_3D_vector_dim2(self):
        """Test 3D vector field with vector dimension 2."""
        vector_field = (
            np.random.rand(self.system3D.xRes, self.system3D.yRes, self.system3D.zRes),
            np.random.rand(self.system3D.xRes, self.system3D.yRes, self.system3D.zRes)
        )
        fig, ax = self.system3D.plot_vector_field(vector_field)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_plot_3D_vector_dim3(self):
        """Test 3D vector field with vector dimension 3."""
        vector_field = (
            np.random.rand(self.system3D.xRes, self.system3D.yRes, self.system3D.zRes),
            np.random.rand(self.system3D.xRes, self.system3D.yRes, self.system3D.zRes),
            np.random.rand(self.system3D.xRes, self.system3D.yRes, self.system3D.zRes)
        )
        fig, ax = self.system3D.plot_vector_field(vector_field)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

if __name__ == '__main__':
    unittest.main()
