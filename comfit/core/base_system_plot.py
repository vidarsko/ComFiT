from typing import Union, Optional

import numpy as np
import scipy as sp

from comfit.tool import tool_colormap
from comfit.tool import tool_create_orthonormal_triad
from comfit.tool import tool_set_plot_axis_properties_matplotlib
from comfit.tool import tool_complete_field
from comfit.tool import tool_print_in_color
from comfit.tool import tool_matplotlib_define_3D_plot_ax
from comfit.tool import tool_plotly_define_3D_plot_ax

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
from comfit.plot import plot_nodes_matplotlib
from comfit.plot import plot_nodes_plotly
from comfit.plot import plot_vector_field_in_plane_both_plot_libs

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

    def plot_prepare(self, field, field_type, **kwargs):

        ## Define figure
        if self.plot_lib == "plotly":
            fig = kwargs.get('fig', go.Figure())
            ax = kwargs.get('ax', {'row': 1, 'col': 1, 'nrows': 1, 'ncols': 1, 'colorbar': False})

        elif self.plot_lib == "matplotlib":
            fig = kwargs['fig'] if 'fig' in kwargs else plt.figure()
            ax = kwargs.get('ax', None)

        # Extend field if not a complete array is given
        if field_type in ['real', 'complex', 'angle']:
            field = tool_complete_field(self, field)
        elif field_type == 'vector':
            field_copy = []
            for n in range(len(field)):
                field_copy.append(tool_complete_field(self, field[n]))
            field = np.array(field_copy)

        ## Checks on field
        # Check if the provided field is bool
        if field.dtype == bool:
            field = field.astype(float)

        # Check if the field is nan
        kwargs['field_is_nan'] = False
        if np.all(np.isnan(field)):
            kwargs['field_is_nan'] = True

        # Check if the field is an order parameter containing several components
        if field_type in ['real', 'complex', 'angle']:
            if field.ndim == self.dim+1:
                tool_print_in_color('Warning: The provided field seems to be an order parameter containing several fields. Only the zeroth component will be plotted.')
                field = field[0]

        # Check if the field is complex when it should be real
        if field_type in ['real', 'vector', 'angle']:
            if np.iscomplexobj(field):
                tool_print_in_color("Warning: the provided field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.")
                print('Max imaginary part: ', np.max(np.imag(field)))
                field = np.real(field)

        ## Define default kwarg values
        kwargs['colormap'] = kwargs.get('colormap', 'viridis')
        kwargs['colormap_object'] = tool_colormap(kwargs['colormap'], plot_lib=self.plot_lib)

        colorbar_default = False
        if field_type == 'complex':
            colorbar_default = True
        if field_type in ['real','angle'] and self.dim>1:
            colorbar_default = True
        if field_type == 'vector' and self.plot_lib == 'plotly':
            colorbar_default = True
        kwargs['colorbar'] = kwargs.get('colorbar',  colorbar_default)

        kwargs['vlim_symmetric'] = kwargs.get('vlim_symmetric', False)

        grid_default = False
        if field_type == 'real' and self.dim in [1,3]:
            grid_default = True
        kwargs['grid'] = kwargs.get('grid', grid_default)

        axis_equal_default = False
        if self.dim in [2,3]:
            axis_equal_default = True
        kwargs['axis_equal'] = kwargs.get('axis_equal', axis_equal_default)

        kwargs['plot_is_3D'] = True if self.dim == 3 else False

        if field_type == 'real':
            # Field limits
            field_min = np.min(field)
            field_max = np.max(field)

            if 'vlim' in kwargs:
                vlim = kwargs['vlim']
                vmin = vlim[0]
                vmax = vlim[1]
            else:
                vmin = field_min
                vmax = field_max
                if 'vlim_symmetric' in kwargs:
                    if kwargs['vlim_symmetric']:
                        vmax = max(abs(vmin), abs(vmax))
                        vmin = -vmax
            
            if self.plot_lib == 'plotly':
                ax['vmin'] = min(vmin, ax.get('vmin', vmin)) 
                ax['vmax'] = max(vmax, ax.get('vmax', vmax))

            else:
                kwargs['vmin'] = vmin
                kwargs['vmax'] = vmax

        return field, fig, ax, kwargs

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
            return plot_complex_field_in_plane_matplotlib(
                self, complex_field, normal_vector, position, **kwargs
            )

    def plot_angle_field_in_plane(
        self,
        angle_field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        **kwargs
    ):
        complex_field = np.exp(1j * angle_field)
        return self.plot_complex_field_in_plane(complex_field, normal_vector=normal_vector, position=position, **kwargs)

    def plot_vector_field_in_plane(
        self,
        vector_field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        spacing = None,
        **kwargs
    ):
        return plot_vector_field_in_plane_both_plot_libs(self, vector_field, normal_vector, position, spacing, **kwargs)

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