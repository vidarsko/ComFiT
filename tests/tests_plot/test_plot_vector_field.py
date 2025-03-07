import unittest
import sys
import os

# Run 
# Adjust the path to import the comfit package
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
import comfit as cf
import numpy as np
import scipy as sp


class TestPlotVectorField(unittest.TestCase):

    def tearDown(self):
        import matplotlib.pyplot as plt
        plt.close('all')

    def test_all_matplotlib(self):
        plot_lib = 'matplotlib'

        bs = cf.BaseSystem(1, plot_lib=plot_lib)
        fig, axs = bs.plot_subplots(3,3)

        for dim in [1,2,3]:
            bs = cf.BaseSystem(dim, plot_lib=plot_lib)

            for n in [1,2,3]:

                if dim == 1 and n == 1:
                    vector_field = [bs.x]
                elif dim == 1 and n == 2:
                    vector_field = [bs.x, 2*bs.x]
                elif dim == 1 and n == 3:
                    vector_field = [bs.x, 2*bs.x, 3*bs.x]
                elif dim == 2 and n == 1:
                    vector_field = [bs.x]
                elif dim == 2 and n == 2:
                    vector_field = [bs.x, bs.y]
                elif dim == 2 and n == 3:
                    vector_field = [bs.x, bs.y, 2*bs.x]
                elif dim == 3 and n == 1:
                    vector_field = [bs.x+bs.y]
                elif dim == 3 and n == 2:
                    vector_field = [bs.x+bs.y, 2*bs.x+2*bs.y]
                elif dim == 3 and n == 3:
                    vector_field = [bs.x+bs.y, 2*bs.x+2*bs.y, 3*bs.x+3*bs.y]

                try:
                    bs.plot_vector_field(vector_field, fig=fig, ax=axs[dim-1][n-1])
                except Exception as e:
                    self.fail(f'plot_vector_field failed for dim={dim}, n={n} with error: {e}')

        # bs.show(fig)

    def test_all_plotly(self):
        plot_lib = 'plotly'

        bs = cf.BaseSystem(1, plot_lib=plot_lib)
        fig, axs = bs.plot_subplots(3,3)

        for dim in [1,2,3]:
            bs = cf.BaseSystem(dim, plot_lib=plot_lib)

            for n in [1,2,3]:

                if dim == 1 and n == 1:
                    vector_field = [bs.x]
                elif dim == 1 and n == 2:
                    vector_field = [bs.x, 2*bs.x]
                elif dim == 1 and n == 3:
                    vector_field = [bs.x, 2*bs.x, 3*bs.x]
                elif dim == 2 and n == 1:
                    vector_field = [bs.x]
                elif dim == 2 and n == 2:
                    vector_field = [bs.x, bs.y]
                elif dim == 2 and n == 3:
                    vector_field = [bs.x, bs.y, 2*bs.x]
                elif dim == 3 and n == 1:
                    vector_field = [bs.x+bs.y]
                elif dim == 3 and n == 2:
                    vector_field = [bs.x+bs.y, 2*bs.x+2*bs.y]
                elif dim == 3 and n == 3:
                    vector_field = [bs.x+bs.y, 2*bs.x+2*bs.y, 3*bs.x+3*bs.y]

                try:
                    bs.plot_vector_field(vector_field, fig=fig, ax=axs[dim-1][n-1])
                except Exception as e:
                    self.fail(f'plot_vector_field failed for dim={dim}, n={n} with error: {e}')

        # bs.show(fig)

    def test_all_plotly_with_colorbar(self):
        plot_lib = 'plotly'

        bs = cf.BaseSystem(1, plot_lib=plot_lib)
        fig, axs = bs.plot_subplots(3,3)

        for dim in [1,2,3]:
            bs = cf.BaseSystem(dim, plot_lib=plot_lib)

            for n in [1,2,3]:

                if dim == 1 and n == 1:
                    vector_field = [bs.x]
                elif dim == 1 and n == 2:
                    vector_field = [bs.x, 2*bs.x]
                elif dim == 1 and n == 3:
                    vector_field = [bs.x, 2*bs.x, 3*bs.x]
                elif dim == 2 and n == 1:
                    vector_field = [bs.x]
                elif dim == 2 and n == 2:
                    vector_field = [bs.x, bs.y]
                elif dim == 2 and n == 3:
                    vector_field = [bs.x, bs.y, 2*bs.x]
                elif dim == 3 and n == 1:
                    vector_field = [bs.x+bs.y]
                elif dim == 3 and n == 2:
                    vector_field = [bs.x+bs.y, 2*bs.x+2*bs.y]
                elif dim == 3 and n == 3:
                    vector_field = [bs.x+bs.y, 2*bs.x+2*bs.y, 3*bs.x+3*bs.y]

                try:
                    bs.plot_vector_field(vector_field, fig=fig, ax=axs[dim-1][n-1], colorbar=True, colormap='sunburst')
                except Exception as e:
                    self.fail(f'plot_vector_field failed for dim={dim}, n={n} with error: {e}')

        # bs.show(fig)

if __name__ == '__main__':
    unittest.main()