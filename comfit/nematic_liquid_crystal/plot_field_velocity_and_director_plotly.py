import numpy as np
import plotly.graph_objects as go
import plotly.figure_factory as ff
from comfit.tool.tool_complete_field import tool_complete_field
from comfit.tool.tool_set_plot_axis_properties_plotly import tool_set_plot_axis_properties_plotly

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

    # Kewyord arguments
    colorbar = kwargs.get('colorbar', True)

    # Extend the field if not a complete array is given
    field = tool_complete_field(self, field)

    fig_tmp = plot_field_plotly(self, field, **kwargs)
    fig.add_trace(fig_tmp.data[0])

    fig_tmp = ff.create_streamline(x=self.x/self.a0, 
                                    y=self.y.flatten()/self.a0, 
                                    u=velocity[0], 
                                    v=velocity[1], 
                                    name='Velocity',
                                    arrow_scale=1)
    fig.add_trace(fig_tmp.data[0])

    kwargs['fig'] = fig
    tool_set_plot_axis_properties_plotly(self, **kwargs)
    return fig