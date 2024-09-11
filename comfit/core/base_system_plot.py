from typing import Union, Optional

import numpy as np
import scipy as sp

from comfit.tool import tool_colormap_angle, tool_colormap_bluewhitered, tool_colormap_sunburst
from comfit.tool import tool_create_orthonormal_triad
from comfit.tool import tool_set_plot_axis_properties_matplotlib
from comfit.tool import tool_complete_field

from comfit.plot import plot_surface_matplotlib
from comfit.plot import plot_field_matplotlib
from comfit.plot import plot_field_plotly
from comfit.plot import plot_complex_field_matplotlib
from comfit.plot import plot_complex_field_plotly
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
    """ Plotting methods for the base system class"""

    def plot_field(self, field, **kwargs):
        return plot_field_plotly(self, field, **kwargs)

    def plot_complex_field(self, complex_field, **kwargs):
        return plot_complex_field_plotly(self, complex_field, **kwargs)

    def plot_angle_field(self, angle_field, **kwargs):
        return plot_angle_field_plotly(self, angle_field, **kwargs)
        
    def plot_vector_field(self, vector_field, **kwargs):
        return plot_vector_field_plotly(self, vector_field, **kwargs)

    def plot_field_in_plane(self, field, **kwargs):
        return plot_field_in_plane_plotly(self, field, **kwargs)
    
    def plot_complex_field_in_plane(self, complex_field, **kwargs):
        return plot_complex_field_in_plane_plotly(self, complex_field, **kwargs)

    def plot_angle_field_in_plane(self, angle_field, **kwargs):
        return plot_angle_field_in_plane_plotly(self, angle_field, **kwargs)

    def plot_vector_field_in_plane(self, vector_field, **kwargs):
        return plot_vector_field_in_plane_plotly(self, vector_field, **kwargs)
