from typing import Union, Dict, Any

# General packages
import plotly.graph_objects as go

# Local packages
from .tool_plotly_find_next_sceneN import tool_plotly_find_next_sceneN

def tool_plotly_define_3D_plot_ax(
        fig: go.Figure, 
        ax: Dict[str, Any]
        ) -> Dict[str, Any]:
    """Defines and validates parameters for a 3D plot axis in Plotly.

    This function checks and sets required parameters for 3D plotting, including
    validating the plot dimension and assigning a scene number.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure object to which the 3D plot will be added.
    ax : dict
        Dictionary containing axis parameters. Must include or will set:
            - plot_dimension : int
                Must be 3 for 3D plots
            - sceneN : int
                Scene number for the 3D plot, auto-generated if not provided
                
    Returns
    -------
    dict
        Updated axis parameters dictionary with validated 3D plot settings.

    Raises
    ------
    ValueError
        If plot_dimension is not 3.
    """

    ax['plot_dimension'] = ax.get('plot_dimension', 3)
    if ax['plot_dimension'] != 3:
        raise ValueError(f"Plot dimension must be 3 for 3D plots. Got {ax['plot_dimension']}")

    ax['sceneN'] = ax.get('sceneN', tool_plotly_find_next_sceneN(fig))

    return ax