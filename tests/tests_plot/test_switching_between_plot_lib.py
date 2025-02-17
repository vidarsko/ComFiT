import unittest
import sys
import os
import numpy as np

import matplotlib
import plotly.graph_objects as go

# Run
# Adjust the path to import the comfit package 
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
import comfit as cf

class TestSwitchingBetweenPlotLib(unittest.TestCase):
    def setUp(self):
        """Initialize the BaseSystem object for testing."""
        self.system1D = cf.BaseSystem(1, xRes=20)
        self.system2D = cf.BaseSystem(2, xRes=20, yRes=20)
        self.system3D = cf.BaseSystem(3, xRes=10, yRes=10, zRes=10)

    def test_plot_field_1d(self):
        # Test 1D field
        field = (np.random.rand(self.system1D.xRes))

        # Single plots
        self.system1D.plot_lib = 'matplotlib'
        fig, ax = self.system1D.plot_field(field)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

        self.system1D.plot_lib = 'plotly'
        fig, ax = self.system1D.plot_field(field)
        self.assertIsInstance(fig, go.Figure)

        # Subplots
        nrows = 2
        ncols = 3

        self.system1D.plot_lib = 'matplotlib'
        fig, axs = self.system1D.plot_subplots(nrows, ncols)
        for i in range(nrows):
            for j in range(ncols):
                fig, ax = self.system1D.plot_field(field, fig=fig, ax=axs[i][j])
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(axs, np.ndarray)

        self.system1D.plot_lib = 'plotly'
        fig, axs = self.system1D.plot_subplots(nrows, ncols)
        for i in range(nrows):
            for j in range(ncols):
                fig, ax = self.system1D.plot_field(field, fig=fig, ax=axs[i][j])
        self.assertIsInstance(fig, go.Figure)

    def test_plot_field_2d(self):
        # Test 2D field
        field = (np.random.rand(self.system2D.xRes, self.system2D.yRes))

        # Single plots
        self.system2D.plot_lib = 'matplotlib'
        fig, ax = self.system2D.plot_field(field)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

        self.system2D.plot_lib = 'plotly'
        fig, ax = self.system2D.plot_field(field)
        self.assertIsInstance(fig, go.Figure)

        # Subplots
        nrows = 2
        ncols = 3

        self.system2D.plot_lib = 'matplotlib'
        fig, axs = self.system2D.plot_subplots(nrows, ncols)
        for i in range(nrows):
            for j in range(ncols):
                fig, ax = self.system2D.plot_field(field, fig=fig, ax=axs[i][j])
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(axs, np.ndarray)

        self.system2D.plot_lib = 'plotly'
        fig, axs = self.system2D.plot_subplots(nrows, ncols)
        for i in range(nrows):
            for j in range(ncols):
                fig, ax = self.system2D.plot_field(field, fig=fig, ax=axs[i][j])
        self.assertIsInstance(fig, go.Figure)


    def test_plot_field_3d(self):
        # Test 3D field
        field = (np.random.rand(self.system3D.xRes, self.system3D.yRes, self.system3D.zRes))

        # Single plots
        self.system3D.plot_lib = 'matplotlib'
        fig, ax = self.system3D.plot_field(field)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

        self.system3D.plot_lib = 'plotly'
        fig, ax = self.system3D.plot_field(field)
        self.assertIsInstance(fig, go.Figure)

        # Subplots
        nrows = 2
        ncols = 3

        self.system3D.plot_lib = 'matplotlib'
        fig, axs = self.system3D.plot_subplots(nrows, ncols)
        for i in range(nrows):
            for j in range(ncols):
                fig, ax = self.system3D.plot_field(field, fig=fig, ax=axs[i][j])
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(axs, np.ndarray)

        self.system3D.plot_lib = 'plotly'
        fig, axs = self.system3D.plot_subplots(nrows, ncols)
        for i in range(nrows):
            for j in range(ncols):
                fig, ax = self.system3D.plot_field(field, fig=fig, ax=axs[i][j])
        self.assertIsInstance(fig, go.Figure)


    
if __name__ == '__main__':
    unittest.main()