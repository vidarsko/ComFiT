import numpy as np
import plotly.graph_objects as go

from typing import Union

from comfit.tool import tool_complete_field
from comfit.tool import tool_set_plot_axis_properties_matplotlib
from comfit.tool import tool_set_plot_axis_properties_plotly
    
from comfit.tool import tool_colormap
from comfit.tool import tool_plotly_find_next_xN
from comfit.tool import tool_plotly_define_2D_plot_ax
from comfit.tool import tool_plotly_define_3D_plot_ax

from comfit.tool import tool_plotly_colorbar

def plot_field_plotly(self, field: np.ndarray, **kwargs) -> go.Figure:
    """Plots the given (real) field.
    
    Args:
        field (array-like): The field to be plotted.
        **kwargs: Keyword arguments for the plot.
            See https://comfitlib.com/ClassBaseSystem/ 
            for a full list of keyword arguments.
    
    Returns:
        The figure containing the plot (plotly.graph_objects.Figure).
    """

    field, fig, ax, kwargs = self.plot_prepare(field, field_type = 'real', **kwargs)

    x = kwargs.get('x', self.x / self.a0)
    y = kwargs.get('y', self.y / self.a0) if self.dim > 1 else None
    z = kwargs.get('z', self.z / self.a0) if self.dim > 2 else None

    ###############################################################
    ###################### DIMENSION: 1 ###########################
    ###############################################################

    if self.dim == 1:

        ax = tool_plotly_define_2D_plot_ax(fig, ax) #Defines xN, yN and plot_dimension

        if not kwargs['field_is_nan']:
            hovertemplate = 'x: %{x:.1f} a₀<br>field: %{y:.1f}'

        if kwargs['fourier']:
            hovertemplate = 'k: %{x:.1f} <br>field: %{y:.1f}'

        if not kwargs['field_is_nan']:
            trace = go.Scatter(
                x=x,
                y=field,
                mode='lines',
                name='',
                hovertemplate=hovertemplate,
                xaxis=ax['xN'],
                yaxis=ax['yN'],
                showlegend=False
            )
            fig.add_trace(trace)

    ###############################################################
    ###################### DIMENSION: 2 ###########################
    ###############################################################

    if self.dim == 2:
        
        ax = tool_plotly_define_2D_plot_ax(fig, ax) #Defines xN, yN and plot_dimension

        X = kwargs.get('X', None)
        Y = kwargs.get('Y', None)

        if X is None or Y is None:
            X, Y = np.meshgrid(x, y, indexing='ij')
            
        if not kwargs['field_is_nan']:
            # Trace
            trace = go.Heatmap(
                x=X.flatten(),
                y=Y.flatten(),
                z=field.flatten(),
                zmin=ax['vmin'],
                zmax=ax['vmax'],
                zsmooth='best',
                hovertemplate='x: %{x:.2f} a₀<br>y: %{y:.2f} a₀<br> field: %{z:.2f}',
                name='',
                colorscale=kwargs['colormap_object'],
                showscale=False,
                xaxis=ax['xN'],
                yaxis=ax['yN']
            )

            fig.add_trace(trace)
        

    ###############################################################
    ###################### DIMENSION: 3 ###########################
    ###############################################################
    elif self.dim == 3:

        ax = tool_plotly_define_3D_plot_ax(fig, ax) #Defines sceneN, plot_dimension

        # Keyword arguments particular to the 3D case
        number_of_layers = kwargs.get('number_of_layers', 1)
        alpha = kwargs.get('alpha', 0.5)

        if 'layer_values' in kwargs:
            layer_values = np.concatenate([[-np.inf], kwargs['layer_values'], [np.inf]])
        else: 
            layer_values = np.linspace(ax['vmin'], ax['vmax'], number_of_layers + 2)


        #Plotting the layers
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

        if not kwargs['field_is_nan']:
            for layer_value in layer_values[1:-1]:
                
                trace = go.Isosurface(
                    x=X.flatten(),
                    y=Y.flatten(),
                    z=Z.flatten(),
                    value = field.flatten(),
                    isomin = layer_value,
                    isomax = layer_value,
                    cmin = ax['vmin'],
                    cmax = ax['vmax'],
                    colorscale = kwargs['colormap_object'],
                    showscale=False,
                    surface=dict(count=3),  # Ensuring only one surface is shown
                    opacity=alpha,
                    scene=ax['sceneN'],
                )

                fig.add_trace(trace)

    if kwargs['colorbar'] and not(ax['colorbar']) and not(kwargs['field_is_nan']):
        ax['colormap_object'] = kwargs['colormap_object']
        fig.add_trace(tool_plotly_colorbar(ax, type='normal'))
        ax['colorbar'] = True

    kwargs['fig'] = fig
    kwargs['ax'] = ax
    tool_set_plot_axis_properties_plotly(self, **kwargs)

    return fig, ax