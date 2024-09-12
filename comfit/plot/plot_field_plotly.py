import numpy as np
import plotly.graph_objects as go

from typing import Union

from comfit.tool import tool_complete_field
from comfit.tool import tool_set_plot_axis_properties_matplotlib
from comfit.tool import tool_set_plot_axis_properties_plotly
    
from comfit.tool import tool_colormap_bluewhitered, tool_colormap_sunburst, tool_complete_field

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

    if field.dtype == bool:
        field = field.astype(float)

    field_is_nan = False
    if np.all(np.isnan(field)):
        field_is_nan = True

    if field.ndim == self.dim+1:
        print("\033[91mWarning - in plot_field: the provided field seems to be an order parameter containing several fields. Only the zeroth component will be plotted.\033[0m")
        field = field[0]

    # Check if the vector field is complex
    if np.iscomplexobj(field):
        print("\033[91mWarning - in plot_field: the provided field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
        print('Max imaginary part: ', np.max(np.imag(field)))
        field = np.real(field)

    fig = kwargs.get('fig', go.Figure())

    # Kewyord arguments
    colorbar = kwargs.get('colorbar', True)
    axis_equal = kwargs.get('axis_equal', None)

    # Extend the field if not a complete array is given
    field = tool_complete_field(self, field)

    kwargs['plot_is_3D'] = False

    # Check if the plot is a subplot
    row = kwargs.get('row', None)
    col = kwargs.get('col', None)

    fig_is_subplot = row is not None and col is not None
    kwargs['fig_is_subplot'] = fig_is_subplot

    ###############################################################
    ###################### DIMENSION: 1 ###########################
    ###############################################################

    if self.dim == 1:

        # Keyword arguments particular to the 1D case
        kwargs['grid'] = kwargs.get('grid', True)

        if axis_equal is None:
            kwargs['axis_equal'] = False

        trace = go.Scatter(
            x=self.x/self.a0,
            y=field,
            mode='lines',
            name='',
            hovertemplate='x: %{x:.2f} a₀<br>field: %{y:.2f}'
        )
        
        if not field_is_nan:
            if fig_is_subplot:
                fig.add_trace(trace, row=row, col=col)
            else:
                fig.add_trace(trace)

    ###############################################################
    ###################### DIMENSION: 2 ###########################
    ###############################################################

    if self.dim == 2:
        
        # Keyword arguments particular to the 2D case
        kwargs['grid'] = kwargs.get('grid', False)


        X, Y = np.meshgrid(self.x, self.y, indexing='ij')
        
        trace = go.Heatmap(
            x=X.flatten()/self.a0,
            y=Y.flatten()/self.a0,
            z=field.flatten(),
            zmin=np.min(field),
            zmax=np.max(field),
            colorscale='Viridis', 
            zsmooth='best',
            hovertemplate='x: %{x:.2f} a₀<br>y: %{y:.2f} a₀<br> field: %{z:.2f}',
            name=''
        )

        if not field_is_nan:
            if fig_is_subplot:
                fig.add_trace(trace, row=row, col=col)
            else:
                fig.add_trace(trace)

    ###############################################################
    ###################### DIMENSION: 3 ###########################
    ###############################################################
    elif self.dim == 3:

        kwargs['plot_is_3D'] = True

        # Keyword arguments particular to the 3D case

        field_min = np.min(field)
        field_max = np.max(field)

        number_of_layers = kwargs.get('number_of_layers', 1)
        alpha = kwargs.get('alpha', 0.5)
        
        if 'vlim' in kwargs:
            vlim = kwargs['vlim']
            vmin = vlim[0]
            vmax = vlim[1]
        else:
            vmin = field_min
            vmax = field_max

        if 'layer_values' in kwargs:
            layer_values = np.concatenate([[-np.inf], kwargs['layer_values'], [np.inf]])
        else: 
            layer_values = np.linspace(vmin, vmax, number_of_layers + 2)

        if 'colorbar' in kwargs:
            colorbar = kwargs['colorbar']
        else:
            colorbar = True

        #Plotting the layers
        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        if not field_is_nan:
            for layer_value in layer_values[1:-1]:
                
                trace = go.Isosurface(
                    x=X.flatten()/self.a0,
                    y=Y.flatten()/self.a0,
                    z=Z.flatten()/self.a0,
                    value = field.flatten(),
                    isomin = layer_value,
                    isomax = layer_value,
                    cmin = vmin,
                    cmax = vmax,
                    surface=dict(count=3),  # Ensuring only one surface is shown
                    colorscale='Viridis',
                    opacity=alpha,
                    showscale=bool(layer_value == layer_values[1])
                )
                if fig_is_subplot:
                    fig.add_trace(trace, row=row, col=col)
                else:
                    fig.add_trace(trace)
    
    kwargs['fig'] = fig
    tool_set_plot_axis_properties_plotly(self,**kwargs)
    return fig