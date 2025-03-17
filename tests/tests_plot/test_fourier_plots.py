import unittest
import sys
import os

# Run 
# Adjust the path to import the comfit package
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
import comfit as cf
import numpy as np
import scipy as sp

show_plots = False

class TestFourierPlots(unittest.TestCase):

    def tearDown(self):
        import matplotlib.pyplot as plt
        plt.close('all')
        
    def test_1D_real_field(self):
        for plot_lib in ['plotly', 'matplotlib']:
            # Initialize BaseSystem object
            bs = cf.BaseSystem(1, plot_lib=plot_lib, xlim=[-10,10])
            # Create field
            field = bs.x/200
            fig, ax = bs.plot_field(field, fourier=True)
            if show_plots:
                bs.show(fig)
            self.assertIsNotNone(fig, msg=f'plot_field with fourier in 1D failed for plot_lib={plot_lib}')
    
    def test_2D_real_field(self):
        for plot_lib in ['plotly', 'matplotlib']:
            # Initialize BaseSystem object
            bs = cf.BaseSystem(2, plot_lib=plot_lib, xlim=[-10,10], ylim=[-10,10])
            # Create field
            field = (bs.x+bs.y)
            fig, ax = bs.plot_field(field, fourier=True)
            if show_plots:
                bs.show(fig)
            self.assertIsNotNone(fig, msg=f'plot_field with fourier in 2D failed for plot_lib={plot_lib}')
    
    def test_3D_real_field(self):
        for plot_lib in ['plotly', 'matplotlib']:
            # Initialize BaseSystem object
            bs = cf.BaseSystem(3, plot_lib=plot_lib, xlim=[-10,10], ylim=[-10,10], zlim=[-10,10])
            # Create field
            field = (bs.x+bs.y+bs.z)
            fig, ax = bs.plot_field(field, fourier=True)
            if show_plots:
                bs.show(fig)
            self.assertIsNotNone(fig, msg=f'plot_field with fourier in 3D failed for plot_lib={plot_lib}')
    
    def test_1D_complex_field(self):
        for plot_lib in ['plotly', 'matplotlib']:
            # Initialize BaseSystem object
            bs = cf.BaseSystem(1, plot_lib=plot_lib, xlim=[-10,10])
            # Create complex field
            field = bs.x*np.exp(1j * bs.x / 5)
            fig, ax = bs.plot_complex_field(field, fourier=True)
            if show_plots:
                bs.show(fig)
            self.assertIsNotNone(fig, msg=f'plot_complex_field with fourier for complex field in 1D failed for plot_lib={plot_lib}')

    def test_2D_complex_field(self):
        for plot_lib in ['plotly', 'matplotlib']:
            # Initialize BaseSystem object
            bs = cf.BaseSystem(2, plot_lib=plot_lib, xlim=[-10,10], ylim=[-10,10])
            # Create complex field
            field = (1 + 0.5*np.sin(bs.x/5)) * np.exp(1j * (bs.x**2 + bs.y**2) / 50)
            fig, ax = bs.plot_complex_field(field, fourier=True)
            if show_plots:
                bs.show(fig)
            self.assertIsNotNone(fig, msg=f'plot_complex_field with fourier for complex field in 2D failed for plot_lib={plot_lib}')

    def test_3D_complex_field(self):
        for plot_lib in ['plotly', 'matplotlib']:
            # Initialize BaseSystem object
            bs = cf.BaseSystem(3, plot_lib=plot_lib, xlim=[-10,10], ylim=[-10,10], zlim=[-10,10])
            # Create complex field
            # Create complex field with amplitude modulation
            field = (bs.x+1j*bs.y+1j*bs.z)/200
            fig, ax = bs.plot_complex_field(field, fourier=True)
            if show_plots:
                bs.show(fig)
            self.assertIsNotNone(fig, msg=f'plot_complex_field with fourier for complex field in 3D failed for plot_lib={plot_lib}')


if __name__ == '__main__':
    unittest.main()