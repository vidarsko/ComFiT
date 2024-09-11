import plotly.graph_objects as go
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from comfit.tool.tool_colormap_angle import tool_colormap_angle


def plot_tool_plotly_add_angle_colorbar3D(self, fig: go.Figure) -> go.Figure:
    """Adds a colorbar to a 3D plot with the angle colormap.

    Args:
        fig: The figure to which the colorbar is added.
    
    Returns:
        The figure with the colorbar added.
    """
    

    return fig