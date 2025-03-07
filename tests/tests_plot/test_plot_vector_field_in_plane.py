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
        fig, axs = bs.plot_subplots(1,3)

        bs = cf.BaseSystem(3, plot_lib=plot_lib)

        for n in [1,2,3]:

            if n == 1:
                vector_field = [bs.x+bs.y]
            elif n == 2:
                vector_field = [bs.x+bs.y, 2*bs.x+2*bs.y]
            elif n == 3:
                vector_field = [bs.x+bs.y, 2*bs.x+2*bs.y, 3*bs.x+3*bs.y]

            try:
                bs.plot_vector_field_in_plane(vector_field, fig=fig, ax=axs[n-1],spacing=1000)
            except Exception as e:
                self.fail(f'plot_vector_field failed for n={n} with error: {e}')

        # bs.show(fig)

    def test_all_plotly(self):
        plot_lib = 'plotly'

        bs = cf.BaseSystem(1, plot_lib=plot_lib)
        fig, axs = bs.plot_subplots(1,3)

        bs = cf.BaseSystem(3, plot_lib=plot_lib)

        for n in [1,2,3]:

            if n == 1:
                vector_field = [bs.x+bs.y]
            elif n == 2:
                vector_field = [bs.x+bs.y, 2*bs.x+2*bs.y]
            elif n == 3:
                vector_field = [bs.x+bs.y, 2*bs.x+2*bs.y, 3*bs.x+3*bs.y]

            try:
                bs.plot_vector_field_in_plane(vector_field, fig=fig, ax=axs[n-1])
            except Exception as e:
                self.fail(f'plot_vector_field_in_plane with plotly failed for n={n} with error: {e}')

        # bs.show(fig)

    # def test_all_plotly_with_colorbar(self):
    #     plot_lib = 'plotly'

    #     bs = cf.BaseSystem(1, plot_lib=plot_lib)
    #     fig, axs = bs.plot_subplots(1,3)


    #     bs = cf.BaseSystem(3, plot_lib=plot_lib)

    #     for n in [1,2,3]:

    #         if n == 1:
    #             vector_field = [bs.x+bs.y]
    #         elif n == 2:
    #             vector_field = [bs.x+bs.y, 2*bs.x+2*bs.y]
    #         elif n == 3:
    #             vector_field = [bs.x+bs.y, 2*bs.x+2*bs.y, 3*bs.x+3*bs.y]

    #         try:
    #             bs.plot_vector_field_in_plane(vector_field, fig=fig, ax=axs[n-1], colorbar=True, colormap='sunburst')
    #         except Exception as e:
    #             self.fail(f'plot_vector_field_in_plane with plotly and custom colorbar failed for n={n} with error: {e}')

        # bs.show(fig)

if __name__ == '__main__':
    unittest.main()