from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

# Standard library imports
from typing import Union

# Third-party library imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go

# Local application imports
from comfit.tool import (
    tool_complete_field,
    tool_set_plot_axis_properties_matplotlib,
    tool_set_plot_axis_properties_plotly,
    tool_colormap,
    tool_matplotlib_define_2D_plot_ax,
    tool_matplotlib_define_3D_plot_ax
)

from .plot_surface_matplotlib import plot_surface_matplotlib

def plot_field_matplotlib(
        self: 'BaseSystem',
        field: np.ndarray,
        **kwargs: Any
        ) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """Plot the given (real) field.

    Parameters
    ----------
    self : BaseSystem
        A BaseSystem (or derived) instance.
    field : np.ndarray
        The field to be plotted.
    kwargs : Any
        Keyword arguments for the plot.
        See https://comfitlib.com/ClassBaseSystem/ 
        for a full list of keyword arguments.

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        - The figure containing the plot
        - The axes containing the plot
    """

    field, fig, ax, kwargs = self.plot_prepare(field, field_type = 'real', **kwargs)

    ###############################################################
    ###################### DIMENSION: 1 ###########################
    ###############################################################

    if self.dim == 1:

        ax = tool_matplotlib_define_2D_plot_ax(fig, ax)

        # Extract coordinates
        x = kwargs.get('x', self.x/self.a0).flatten()

        ax.plot(self.x/self.a0, field)

    ###############################################################
    ###################### DIMENSION: 2 ###########################
    ###############################################################

    if self.dim == 2:
        
        # Keyword arguments particular to the 2D case
        kwargs['grid'] = kwargs.get('grid', False)

        ax = tool_matplotlib_define_2D_plot_ax(fig, ax)
            
        # Set the colormap
        colormap_string = kwargs.get('colormap', 'viridis')
        colormap = tool_colormap(colormap_string, plot_lib='matplotlib')
        
        # Value limits symmetric
        vlim_symmetric = kwargs.get('vlim_symmetric', False)

        # Extract coordinates
        x = kwargs.get('x', self.x/self.a0).flatten()
        dx = x[1] - x[0]
        xmin = x[0]
        xmax = x[-1]+dx

        y = kwargs.get('y', self.y/self.a0).flatten()
        dy = y[1] - y[0]
        ymin = y[0]
        ymax = y[-1]+dy

        X = kwargs.get('X', None)
        Y = kwargs.get('Y', None)

        if X is None or Y is None:
            X, Y = np.meshgrid(x, y, indexing='ij')

        pcm = ax.pcolormesh(X , Y , field, shading='gouraud', cmap=colormap)

        xlim = [xmin, xmax-dx]
        ylim = [ymin, ymax-dy]

        limits_provided = False
        if 'xlim' in kwargs:
            xlim = np.array(kwargs['xlim'])/self.a0
            limits_provided = True
        else:
            if 'xmin' in kwargs:
                xlim[0] = kwargs['xmin']/self.a0
                limits_provided = True
            
            if 'xmax' in kwargs:
                xlim[1] = kwargs['xmax']/self.a0
                limits_provided = True

        if 'ylim' in kwargs:
            ylim = np.array(kwargs['ylim'])/self.a0
            limits_provided = True
        else:
            if 'ymin' in kwargs:
                ylim[0] = kwargs['ymin']/self.a0
                limits_provided = True
                
            if 'ymax' in kwargs:
                ylim[1] = kwargs['ymax']/self.a0
                limits_provided = True

        # If explicit limits are provided, use them to change the vlim ranges
        if limits_provided:
            region_to_plot = np.zeros(self.dims).astype(bool)
            region_to_plot[(xlim[0] <= X)*(X <= xlim[1])*(ylim[0] <= Y)*(Y <= ylim[1])] = True
            vlim = [np.min(field[region_to_plot]), np.max(field[region_to_plot])]

        else:
            vlim = [np.min(field), np.max(field)]
        
        # Set the value limits
        if 'vlim' in kwargs:
            vlim = kwargs['vlim']
        else:
            if 'vmin' in kwargs:
                vlim[0] = kwargs['vmin']
            if 'vmax' in kwargs:
                vlim[1] = kwargs['vmax']

        if vlim[1] - vlim[0] < 1e-10:
            vlim = [vlim[0]-0.05, vlim[1]+0.05]

        pcm.set_clim(vmin=vlim[0], vmax=vlim[1])

        if 'vlim_symmetric' in kwargs:
            vlim_symmetric = kwargs['vlim_symmetric']
            if vlim_symmetric:
                cmax = abs(field).max()
                cmin = -cmax
                pcm.set_clim(vmin=cmin, vmax=cmax)

        colorbar = kwargs.get('colorbar', True)

        if colorbar:
            cbar = plt.colorbar(pcm, ax=ax)

    ###############################################################
    ###################### DIMENSION: 3 ###########################
    ###############################################################

    elif self.dim == 3:

        ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
        kwargs['plot_is_3D'] = True
            
        # Keyword arguments particular to the 3D case
        number_of_layers = kwargs.get('number_of_layers', 1)
        alpha = kwargs.get('alpha', 0.5)

        vmin = kwargs['vmin']
        vmax = kwargs['vmax']

        if 'layer_values' in kwargs:
            layer_values = np.concatenate([[-np.inf], kwargs['layer_values'], [np.inf]])
        else: 
            layer_values = np.linspace(vmin, vmax, number_of_layers + 2)

        field_min = np.min(field)
        field_max = np.max(field)

        if field_min < layer_values[1] < field_max:
            plot_surface_matplotlib(self, field=field, 
                                    value=layer_values[1], 
                                    color=kwargs['colormap_object']((layer_values[1]-vmin) / (vmax-vmin)), 
                                    alpha=alpha,
                                    ax=ax)

        for layer_value in layer_values[2:-1]:
            if field_min < layer_value < field_max:
                plot_surface_matplotlib(self, field=field, 
                                        value=layer_value, 
                                        color=kwargs['colormap_object']((layer_value-vmin) / (vmax-vmin)), 
                                        alpha=alpha,
                                        ax=ax)

        if kwargs['colorbar']:
            sm = plt.cm.ScalarMappable(cmap=kwargs['colormap_object'])
            sm.set_clim(kwargs['vmin'], kwargs['vmax'])
            plt.colorbar(sm, ax=ax, pad=0.2)

    kwargs['ax'] = ax
    ax = tool_set_plot_axis_properties_matplotlib(self, **kwargs)
    return fig, ax
