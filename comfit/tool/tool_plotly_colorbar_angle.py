import numpy as np
import plotly.graph_objects as go
from comfit.tool import tool_colormap_angle


def tool_plotly_colorbar_angle(ax=None):
    """Add a colorbar for the angle field to a plotly figure.

    Args:
        ax: The (row,column) index of the subplot to add the colorbar to.
    
    Returns:
        plotly.graph_objects.Scatter: A scatter plot with no data that can be used as a colorbar.
    """

    custom_colormap = tool_colormap_angle(pyplot=True)

    trace = go.Scatter(
            x=[None], y=[None], mode='markers',
            showlegend=False,
            marker=dict(
                colorscale=custom_colormap,
                cmin=-np.pi, cmax=np.pi,
                colorbar=dict(
                tickvals=[-np.pi, -2*np.pi/3, -np.pi/3, 0, np.pi/3, 2*np.pi/3, np.pi],
                ticktext=['-π', '-2π/3', '-π/3', '0', 'π/3', '2π/3', 'π']
                )
            ),
            hoverinfo='none'
            )
    
    if ax is not None:
        # The idea here is to place the colorbar in the middle of the subplot
        # It is not perfect yet (Vidar, 13.02.25)
        trace.marker.colorbar.update(
                len=1/ax[0,1],
                xanchor='right',
                x=ax[1,0]/ax[1,1],
                yanchor='middle',
                y=1-1/ax[0,1]/2-(ax[0,0]-1)/ax[0,1]
            )

    return trace