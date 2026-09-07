from typing import Dict, Any

# General packages
import plotly.graph_objects as go

# Local packages
from comfit.tool import tool_plotly_find_next_xN

def tool_plotly_define_2D_plot_ax(
        fig: go.Figure, 
        ax: Dict[str, Any]
        ) -> Dict[str, Any]:
    """Define x-axis, y-axis and plot dimension parameters for 2D Plotly plots.

    This function sets up axis naming conventions and validates plot dimensions
    for 2D plots in Plotly. It determines the next available x-axis identifier
    and corresponding y-axis.
    
    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure object to define axes for.
    ax : Dict[str, Any]
        Dictionary containing axis parameters. May include:
        xN : str (optional, X-axis identifier, e.g., 'x', 'x2', 'x3'),
        yN : str (optional, Y-axis identifier, e.g., 'y', 'y2', 'y3'),
        plot_dimension : int (optional, plot dimension, must be 2 for 2D plots).
    
    Returns
    -------
    Dict[str, Any]
        Updated dictionary containing:
        xN : str (X-axis identifier),
        yN : str (Y-axis identifier),
        plot_dimension : int (Plot dimension, always 2),
        xaxisN : str (Full x-axis reference, e.g., 'xaxis', 'xaxis2'),
        yaxisN : str (Full y-axis reference, e.g. 'yaxis', 'yaxis2').

    Raises
    ------
    ValueError
        If plot_dimension is not 2.

    Notes
    -----
    The function automatically determines the next available x-axis identifier
    if not provided in the input dictionary.
    """

    #Plot nature
    ax['xN'] = ax.get('xN', tool_plotly_find_next_xN(fig))
    ax['yN'] = 'y' if ax['xN'] == 'x' else f'y{ax["xN"][1:]}'

    # Plot dimension
    ax['plot_dimension'] = ax.get('plot_dimension', 2)
    if ax['plot_dimension'] != 2:
        raise ValueError(f"Plot dimension must be 2 for 2D plots. Got {ax['plot_dimension']}")

    # Axis number
    axis_number = 1 if ax['xN'] == 'x' else int(ax['xN'][1:])  # extracts everything after 'x' and converts to integer
    
    ax['xaxisN'] = f'xaxis{axis_number}'
    ax['yaxisN'] = f'yaxis{axis_number}'

    return ax