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

from comfit.plot import plot_save_matplotlib
from comfit.plot import plot_save_plotly

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

    def plot_nodes(self, disclination_nodes, **kwargs):
        if self.plot_lib == 'plotly':
            return plot_nodes_plotly(self, disclination_nodes, **kwargs)
        elif self.plot_lib == 'matplotlib':
            return plot_nodes_matplotlib(self, disclination_nodes, **kwargs)


    # Figure handling methods
    def plot_subplots(self, number_of_rows, number_of_columns, *args, **kwargs):
        if self.plot_lib == "plotly":
            return plot_subplots_plotly(number_of_rows, number_of_columns, *args, **kwargs)
        elif self.plot_lib == "matplotlib":
            return plot_subplots_matplotlib(number_of_rows, number_of_columns, *args, **kwargs)

    def plot_save(self, counter, fig, **kwargs):
        if self.plot_lib == "plotly":
            return plot_save_plotly(counter, fig, **kwargs)
        elif self.plot_lib == "matplotlib":
            return plot_save_matplotlib(counter, fig, **kwargs)

    def show(self, fig):
        if self.plot_lib == "matplotlib":
            plt.show()
        elif self.plot_lib == "plotly":
            fig.show()