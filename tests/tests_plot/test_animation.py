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

    def test_save_plotly(self):
        bs = cf.BaseSystem(dim=1, xmin=0, xmax=10, xRes=100, plot_lib='plotly')

        fig, ax = bs.plot_field(np.sin(bs.x))
        
        try:
            bs.plot_save(fig)
        except:
            self.fail('plot_save_plotly raised an exception.')
        
    def test_save_matplotlib(self):
        bs = cf.BaseSystem(dim=1, xmin=0, xmax=10, xRes=100, plot_lib='matplotlib')

        fig, ax = bs.plot_field(np.sin(bs.x))
        
        try:
            bs.plot_save(fig)
        except:
            self.fail('plot_save_matplotlib raised an exception.')


    def test_animation_simple_plotly(self):
        
        bs = cf.BaseSystem(dim=1, xmin=0, xmax=10, xRes=100, plot_lib='plotly')
        
        for n in range(20):
            fig, ax = bs.plot_field(np.sin(bs.x-n/20))
            bs.plot_save(fig, n)
        
        try:
            cf.tool_make_animation_movie(n, delete_png=False)
        except:
            self.fail('tool_make_animation_movie_plotly raised an exception.')

        try:
            cf.tool_make_animation_gif(n, delete_png=True)
        except:
            self.fail('tool_make_animation_gif_plotly raised an exception.')

    
    def test_animation_simple_matplotlib(self):
        
        bs = cf.BaseSystem(dim=1, xmin=0, xmax=10, xRes=100, plot_lib='matplotlib')
        
        for n in range(20):
            fig, ax = bs.plot_field(np.sin(bs.x-n/20))
            bs.plot_save(fig, n)
        
        try:
            cf.tool_make_animation_movie(n, delete_png=False)
        except:
            self.fail('tool_make_animation_movie_matplotlib raised an exception.')

        try:
            cf.tool_make_animation_gif(n, delete_png=True)
        except:
            self.fail('tool_make_animation_gif_matplotlib raised an exception.')
    

    def test_animation_with_ID_plotly(self):

        bs = cf.BaseSystem(dim=1, xmin=0, xmax=10, xRes=100)

        for n in range(20):
            fig, ax = bs.plot_field(np.sin(bs.x-n/20))
            bs.plot_save(fig, n, ID='test')
        
        try:
            cf.tool_make_animation_movie(n, ID='test', delete_png=False)
        except:
            self.fail('test_animation_with_ID_plotly raised an exception.')

        try:
            cf.tool_make_animation_gif(n, ID='test', delete_png=True)
        except:
            self.fail('test_animation_with_ID_plotly raised an exception.')
        
    def test_animation_with_ID_matplotlib(self):

        bs = cf.BaseSystem(dim=1, xmin=0, xmax=10, xRes=100, plot_lib='matplotlib')

        for n in range(20):
            fig, ax = bs.plot_field(np.sin(bs.x-n/20))
            bs.plot_save(fig, n, ID='test')

        try:
            cf.tool_make_animation_movie(n, ID='test', delete_png=False)
        
        except:
            self.fail('test_animation_with_ID_matplotlib raised an exception.')




if __name__ == '__main__':
    unittest.main()