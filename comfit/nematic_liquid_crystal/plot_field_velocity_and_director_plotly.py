import numpy as np
import plotly.graph_objects as go
import plotly.figure_factory as ff
from comfit.tool.tool_complete_field import tool_complete_field
from comfit.tool.tool_set_plot_axis_properties_plotly import tool_set_plot_axis_properties_plotly
from comfit.tool.tool_add_spacing_2D import tool_add_spacing_2D

from comfit.plot.plot_field_plotly import plot_field_plotly


def plot_field_velocity_and_director_plotly(self, field, velocity, director, **kwargs):
    """Plot the fields, velocity, and director field in 2 dimensions using Plotly.

    Args:
        field (ndarray): The field to be plotted.
        velocity (ndarray): The velocity to be plotted.
        director (ndarray): The director to be plotted.
        **kwargs: Keyword arguments for the plot.

    Returns:
        The plotly figure (go.Figure).
    """
    if field.dtype == bool:
        field = field.astype(float)

    # Check if the vector field is complex
    if np.iscomplexobj(field):
        print(
                "\033[91mWarning: the provided field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
        print('Max imaginary part: ', np.max(np.imag(field)))
    field = np.real(field)

    # Check if an axis object is provided
    fig = kwargs.get('fig', go.Figure())

    # Extend the field if not a complete array is given
    field = tool_complete_field(self, field)

    kwargs['colormap'] = kwargs.get('colormap', 'Picnic')
    kwargs['vlim_symmetric'] = kwargs.get('vlim_symmetric', True)

    fig_tmp = plot_field_plotly(self, field, **kwargs)
    fig.add_trace(fig_tmp.data[0].update(showlegend=False))

    # Plot the director field using quiver
    X, Y = np.meshgrid(self.x, self.y, indexing='ij')

    spacing = kwargs.get('spacing', 5)

    X, Y, U, V = tool_add_spacing_2D(X,Y,director[0],director[1],spacing)

    u = U.flatten()
    v = V.flatten()

    magnitude = np.sqrt(u**2 + v**2)
    magnitude_max = max(np.max(magnitude),1e-12)
    magnitude_normalized = magnitude/magnitude_max

    angle = np.arctan2(v, u)
    direction = np.array([np.cos(angle), np.sin(angle)]).T
    
    fig.add_trace(go.Scatter(
        x=X.flatten()/self.a0,
        y=Y.flatten()/self.a0,
        mode='markers',
        marker=dict(symbol='line-ew', 
            angle=90-angle.flatten()*180/np.pi, 
            size=2*spacing*magnitude_normalized.flatten(), 
            sizemode='diameter',
            color=magnitude.flatten(), 
            cmin=0,
            cmax=magnitude_max,
            line=dict(color='black')
            ),
            hovertemplate='<b>x:</b> %{x:.2f}a0<br>' +
                '<b>y:</b> %{y:.2f}a0<br>' +
                '<b>nx:</b> %{customdata[0]:.2e}<br>' +  
                '<b>ny:</b> %{customdata[1]:.2e}<br>',
            customdata=np.stack((u.flatten(), v.flatten()), axis=-1), showlegend=False  # Adding ux, uy and u as customdata
        )
        )

    fig_tmp = ff.create_streamline(x=self.x/self.a0, 
                                    y=self.y.flatten()/self.a0, 
                                    u=velocity[0].T, 
                                    v=velocity[1].T, 
                                    arrow_scale=1, 
                                    density=1)
        
    for trace in fig_tmp.data:
        trace.showlegend = False

    # fig.add_trace(fig_tmp.data[0])

    kwargs['fig'] = fig
    tool_set_plot_axis_properties_plotly(self, **kwargs)
    return fig