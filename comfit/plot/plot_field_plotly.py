import numpy as np
import plotly.graph_objects as go

from typing import Union

from comfit.tool import tool_complete_field
from comfit.tool import tool_set_plot_axis_properties_matplotlib
from comfit.tool import tool_set_plot_axis_properties_plotly
    
from comfit.tool import tool_colormap_bluewhitered, tool_colormap_sunburst, tool_complete_field
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
    ax = kwargs.get('ax', {'row': 1, 'col': 1, 'nrows': 1, 'ncols': 1})

    # Kewyord arguments
    axis_equal = kwargs.get('axis_equal', None)

    # Extend the field if not a complete array is given
    field = tool_complete_field(self, field)

    kwargs['plot_is_3D'] = False

    # Colormap
    colormap = kwargs.get('colormap', 'Viridis')

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
    
    ax['vmin'] = min(vmin, ax.get('vmin', vmin)) 
    ax['vmax'] = max(vmax, ax.get('vmax', vmax))

    ###############################################################
    ###################### DIMENSION: 1 ###########################
    ###############################################################

    if self.dim == 1:

        ax = tool_plotly_define_2D_plot_ax(ax, fig) #Defines xN, yN and plot_dimension

        # Keyword arguments particular to the 1D case
        kwargs['grid'] = kwargs.get('grid', True)
        kwargs['colorbar'] = kwargs.get('colorbar', False)

        if axis_equal is None:
            kwargs['axis_equal'] = False

        if not field_is_nan:
            trace = go.Scatter(
                x=self.x/self.a0,
                y=field,
                mode='lines',
                name='',
                hovertemplate='x: %{x:.1f} a₀<br>field: %{y:.1f}',
                xaxis=ax['xN'],
                yaxis=ax['yN'],
                showlegend=False
            )
            fig.add_trace(trace)

    ###############################################################
    ###################### DIMENSION: 2 ###########################
    ###############################################################

    if self.dim == 2:
        
        ax = tool_plotly_define_2D_plot_ax(ax, fig) #Defines xN, yN and plot_dimension

        # Plot specific keyword arguments
        kwargs['grid'] = kwargs.get('grid', False)
        kwargs['colorbar'] = kwargs.get('colorbar', True)


        X = kwargs.get('X', None)
        Y = kwargs.get('Y', None)

        if X is None or Y is None:
            X, Y = np.meshgrid(self.x, self.y, indexing='ij')
            
        if not field_is_nan:
            # Trace
            trace = go.Heatmap(
                x=X.flatten()/self.a0,
                y=Y.flatten()/self.a0,
                z=field.flatten(),
                zmin=vmin,
                zmax=vmax,
                zsmooth='best',
                hovertemplate='x: %{x:.2f} a₀<br>y: %{y:.2f} a₀<br> field: %{z:.2f}',
                name='',
                colorscale=colormap,
                showscale=False,
                xaxis=ax['xN'],
                yaxis=ax['yN']
            )

            fig.add_trace(trace)
        

    ###############################################################
    ###################### DIMENSION: 3 ###########################
    ###############################################################
    elif self.dim == 3:

        ax = tool_plotly_define_3D_plot_ax(ax, fig) #Defines sceneN, plot_dimension


        # Keyword arguments particular to the 3D case
        number_of_layers = kwargs.get('number_of_layers', 1)
        alpha = kwargs.get('alpha', 0.5)
        kwargs['colorbar'] = kwargs.get('colorbar', True)

        if 'layer_values' in kwargs:
            layer_values = np.concatenate([[-np.inf], kwargs['layer_values'], [np.inf]])
        else: 
            layer_values = np.linspace(vmin, vmax, number_of_layers + 2)


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
                    colorscale = colormap,
                    showscale=False,
                    surface=dict(count=3),  # Ensuring only one surface is shown
                    opacity=alpha,
                    scene=ax['sceneN'],
                )

                fig.add_trace(trace)
    
    if kwargs['colorbar'] and not field_is_nan:
        fig.add_trace(tool_plotly_colorbar(ax, type='normal'))

    kwargs['fig'] = fig
    kwargs['ax'] = ax
    tool_set_plot_axis_properties_plotly(self, **kwargs)
    return fig, ax