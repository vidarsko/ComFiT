from typing import Union, Optional

import numpy as np
import scipy as sp

from comfit.tool import (
    tool_colormap_angle,
    tool_colormap_bluewhitered,
    tool_colormap_sunburst,
)
from comfit.tool import tool_create_orthonormal_triad
from comfit.tool import tool_set_plot_axis_properties_matplotlib
from comfit.tool import tool_complete_field
from comfit.tool import tool_print_in_color

from comfit.plot import plot_surface_matplotlib
from comfit.plot import plot_field_matplotlib
from comfit.plot import plot_field_plotly
from comfit.plot import plot_complex_field_matplotlib
from comfit.plot import plot_complex_field_plotly
from comfit.plot import plot_angle_field_matplotlib
from comfit.plot import plot_angle_field_plotly
from comfit.plot import plot_vector_field_matplotlib
from comfit.plot import plot_vector_field_plotly
from comfit.plot import plot_field_in_plane_matplotlib
from comfit.plot import plot_field_in_plane_plotly
from comfit.plot import plot_complex_field_in_plane_matplotlib
from comfit.plot import plot_complex_field_in_plane_plotly
from comfit.plot import plot_angle_field_in_plane_matplotlib
from comfit.plot import plot_angle_field_in_plane_plotly
from comfit.plot import plot_vector_field_in_plane_matplotlib
from comfit.plot import plot_vector_field_in_plane_plotly
from comfit.plot import plot_nodes_matplotlib
from comfit.plot import plot_nodes_plotly

from comfit.plot import plot_subplots_matplotlib
from comfit.plot import plot_subplots_plotly

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.tri as mtri
import matplotlib.cm as cm
import matplotlib.axes
import matplotlib.figure

import mpl_toolkits

import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.colors as pc

from skimage.measure import marching_cubes

import time


class BaseSystemPlot:
    """Plotting methods for the base system class"""

    def plot_field(self, field: np.ndarray, **kwargs):
        if self.plot_lib == "plotly":
            return plot_field_plotly(self, field, **kwargs)
        elif self.plot_lib == "matplotlib":
            return plot_field_matplotlib(self, field, **kwargs)

    def plot_complex_field(self, complex_field: np.ndarray, **kwargs):
        if self.plot_lib == "plotly":
            return plot_complex_field_plotly(self, complex_field, **kwargs)
        elif self.plot_lib == "matplotlib":
            return plot_complex_field_matplotlib(self, complex_field, **kwargs)

    def plot_angle_field(self, angle_field: np.ndarray, **kwargs):
        if self.plot_lib == "plotly":
            return plot_angle_field_plotly(self, angle_field, **kwargs)
        elif self.plot_lib == "matplotlib":
            return plot_angle_field_matplotlib(self, angle_field, **kwargs)

    def plot_vector_field(self, vector_field: np.ndarray, **kwargs):
        if self.plot_lib == "plotly":
            return plot_vector_field_plotly(self, vector_field, **kwargs)
        elif self.plot_lib == "matplotlib":
            return plot_vector_field_matplotlib(self, vector_field, **kwargs)

    def plot_field_in_plane(
        self,
        field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        **kwargs
    ):
        if self.plot_lib == "plotly":
            return plot_field_in_plane_plotly(
                self, field, normal_vector, position, **kwargs
            )
        elif self.plot_lib == "matplotlib":
            return plot_field_in_plane_matplotlib(
                self, field, normal_vector, position, **kwargs
            )

    def plot_complex_field_in_plane(
        self,
        complex_field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        **kwargs
    ):
        if self.plot_lib == "plotly":
            return plot_complex_field_in_plane_plotly(
                self, complex_field, normal_vector, position, **kwargs
            )
        elif self.plot_lib == "matplotlib":
            return plot_complex_field_in_plane_plotly(
                self, complex_field, normal_vector, position, **kwargs
            )

    def plot_angle_field_in_plane(
        self,
        angle_field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        **kwargs
    ):
        if self.plot_lib == "plotly":
            return plot_angle_field_in_plane_plotly(
                self, angle_field, normal_vector, position, **kwargs
            )
        elif self.plot_lib == "matplotlib":
            return plot_angle_field_in_plane_plotly(
                self, angle_field, normal_vector, position, **kwargs
            )

    def plot_vector_field_in_plane(
        self,
        vector_field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        spacing: int = 2,
        **kwargs
    ):
        if self.plot_lib == "plotly":
            return plot_vector_field_in_plane_plotly(
                self, vector_field, normal_vector, position, spacing, **kwargs
            )
        elif self.plot_lib == "matplotlib":
            return plot_vector_field_in_plane_matplotlib(
                self, vector_field, normal_vector, position, spacing, **kwargs
            )

    def plot_nodes(self, nodes, **kwargs):
        if self.plot_lib == 'plotly':
            return plot_nodes_plotly(self, nodes, **kwargs)
        elif self.plot_lib == 'matplotlib':
            return plot_nodes_matplotlib(self, nodes, **kwargs)


    # Figure handling methods
    def plot_subplots(self, number_of_rows, number_of_columns, **kwargs):
        if self.plot_lib == "plotly":
            return plot_subplots_plotly(number_of_rows, number_of_columns, **kwargs)
        elif self.plot_lib == "matplotlib":
            return plot_subplots_matplotlib(number_of_rows, number_of_columns, **kwargs)

    def plot_save(self, fig, counter=None, **kwargs):

        # Check if the figure is the first argument (deprecated)
        if isinstance(fig, int):
            tool_print_in_color("Warning: The plot_save method has been called with the figure as the second argument. This is deprecated. Please call the method with the figure as the first argument, then counter.")
            counter_tmp = fig
            fig = counter
            counter = counter_tmp
        
        # Get ID
        ID=kwargs.get('ID', None)

        # Keyword arguments
        image_size_inches=kwargs.get('image_size_inches', (6,5))
        dpi=kwargs.get('dpi', 100)

        # Set filename
        if counter is None:
            if ID is None:
                filename = 'plot.png'
            else:
                filename = f'plot_{ID}.png'
        else:
            if ID is None:
                filename = f'plot_{counter}.png'
            else:
                filename = f'plot_{counter}_{ID}.png'


        # Save the figure
        if self.plot_lib == "plotly":
            fig.write_image(filename)

        elif self.plot_lib == "matplotlib":
            fig.set_size_inches(image_size_inches)
            fig.savefig(filename, dpi=dpi)


    def show(self, fig):
        if self.plot_lib == "matplotlib":
            plt.show()
        elif self.plot_lib == "plotly":
            fig.show()