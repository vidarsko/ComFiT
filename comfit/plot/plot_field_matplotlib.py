import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

from typing import Union

from comfit.tool import tool_complete_field
from comfit.tool import tool_set_plot_axis_properties_matplotlib
from comfit.tool import tool_set_plot_axis_properties_plotly

from comfit.tool import tool_colormap
from comfit.tool import tool_matplotlib_define_2D_plot_ax
from comfit.tool import tool_matplotlib_define_3D_plot_ax

from comfit.plot.plot_surface_matplotlib import plot_surface_matplotlib

import matplotlib
from mpl_toolkits.mplot3d import Axes3D

def plot_field_matplotlib(self, 
        field: np.ndarray, 
        **kwargs
        ) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """Plots the given (real) field.
    
    Args:
        field (array-like): The field to be plotted.
        **kwargs: Keyword arguments for the plot.
            See https://comfitlib.com/ClassBaseSystem/ 
            for a full list of keyword arguments.
    
    Returns:
        tuple: A tuple containing:
            - matplotlib.figure.Figure: The figure containing the plot.
            - matplotlib.axes.Axes: The axes containing the plot
    """

    field, fig, ax, kwargs = self.plot_prepare(field, field_type = 'real', **kwargs)

    ###############################################################
    ###################### DIMENSION: 1 ###########################
    ###############################################################

    if self.dim == 1:

        ax = tool_matplotlib_define_2D_plot_ax(fig, ax)

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

        X = kwargs.get('X', None)
        Y = kwargs.get('Y', None)

        if X is None or Y is None:
            X, Y = np.meshgrid(self.x, self.y, indexing='ij')

        pcm = ax.pcolormesh(X / self.a0, Y / self.a0, field, shading='gouraud', cmap=colormap)

        xlim = [self.xmin, self.xmax-self.dx]
        ylim = [self.ymin, self.ymax-self.dy]

        limits_provided = False
        if 'xlim' in kwargs:
            xlim = kwargs['xlim']
            limits_provided = True
        else:
            if 'xmin' in kwargs:
                xlim[0] = kwargs['xmin']
                limits_provided = True
            
            if 'xmax' in kwargs:
                xlim[1] = kwargs['xmax']
                limits_provided = True

        if 'ylim' in kwargs:
            ylim = kwargs['ylim']
            limits_provided = True
        else:
            if 'ymin' in kwargs:
                ylim[0] = kwargs['ymin']
                limits_provided = True
                
            if 'ymax' in kwargs:
                ylim[1] = kwargs['ymax']
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
