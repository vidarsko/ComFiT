import unittest
import sys
import os

# Run 
# Adjust the path to import the comfit package
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
import comfit as cf
import numpy as np
import scipy as sp


class TestPlotAngleField(unittest.TestCase):

  def tearDown(self):
      import matplotlib.pyplot as plt
      plt.close('all')

  def test_3D(self):
    for plot_lib in ['plotly', 'matplotlib']:
      # Initialize BaseSystem object
      bs = cf.BaseSystem(3, plot_lib=plot_lib, xlim=[-10,10], ylim=[-10,10], zlim=[-10,10])
      # Create field
      field = (bs.x+bs.y+bs.z)
      fig, ax = bs.plot_angle_field_in_plane(field)
      # bs.show(fig)
      self.assertIsNotNone(fig, msg=f'plot_angle_field_in_plane in 3D failed for plot_lib={plot_lib}')

  def test_colorbar_3D(self):
    for plot_lib in ['plotly', 'matplotlib']:
      # Initialize BaseSystem object
      bs = cf.BaseSystem(3, plot_lib=plot_lib, xlim=[-10,10], ylim=[-10,10], zlim=[-10,10])
      colormap = 'sunburst'
      # Create field
      field = (bs.x+bs.y+bs.z)
      fig, ax = bs.plot_angle_field_in_plane(field, colorbar=True, colormap=colormap)
      # bs.show(fig)
      self.assertIsNotNone(fig, msg=f'plot_angle_field_in_plane in 3D with custom colorbar failed for plot_lib={plot_lib}')

if __name__ == '__main__':
  unittest.main()