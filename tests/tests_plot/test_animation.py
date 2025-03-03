import unittest
import sys
import os

# Run 
# Adjust the path to import the comfit package
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
import comfit as cf
import numpy as np
import scipy as sp

class TestAnimation(unittest.TestCase):
    def test_animation_simple(self):
        
        bs = cf.BaseSystem(dim=1, xmin=0, xmax=10, xRes=100)
        
        for n in range(20):
            fig, ax = bs.plot_field(np.sin(bs.x-n/20))
            bs.plot_save(n, fig)
        
        try:
            cf.tool_make_animation_movie(n, delete_png=False)
        except:
            self.fail('tool_make_animation_movie raised an exception.')

        try:
            cf.tool_make_animation_gif(n, delete_png=False)
        except:
            self.fail('tool_make_animation_gif raised an exception.')
    




if __name__ == '__main__':
    unittest.main()