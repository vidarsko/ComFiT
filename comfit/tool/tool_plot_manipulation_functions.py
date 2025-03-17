from typing import Optional

# General packages
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def tool_zoom_plot_matplotlib(
        zoom_factor: float, 
        ax: Optional[plt.Axes] = None
        ) -> None:
    """Zoom in/out on a matplotlib plot by adjusting axis limits.

    Parameters
    ----------
    zoom_factor : float
        Factor by which to zoom. Values > 1 zoom in, values < 1 zoom out.
    ax : plt.Axes, optional
        The axes to apply zoom to. If None, uses current axes, by default None.

    Returns
    -------
    None
        Modifies the plot in place.
    """
    if ax is None:
        ax = plt.gca()

    if isinstance(ax, Axes3D):
        # Zoom for 3D plot
        x_lim = ax.get_xlim()
        y_lim = ax.get_ylim()
        z_lim = ax.get_zlim()

        x_center = (x_lim[0] + x_lim[1]) / 2
        y_center = (y_lim[0] + y_lim[1]) / 2
        z_center = (z_lim[0] + z_lim[1]) / 2

        ax.set_xlim3d([x_center - (x_center - x_lim[0]) / zoom_factor,
                       x_center + (x_lim[1] - x_center) / zoom_factor])
        ax.set_ylim3d([y_center - (y_center - y_lim[0]) / zoom_factor,
                       y_center + (y_lim[1] - y_center) / zoom_factor])
        ax.set_zlim3d([z_center - (z_center - z_lim[0]) / zoom_factor,
                       z_center + (z_lim[1] - z_center) / zoom_factor])
    else:
        # Zoom for 2D plot
        x_lim = ax.get_xlim()
        y_lim = ax.get_ylim()

        x_center = (x_lim[0] + x_lim[1]) / 2
        y_center = (y_lim[0] + y_lim[1]) / 2

        ax.set_xlim([x_center - (x_center - x_lim[0]) / zoom_factor,
                     x_center + (x_lim[1] - x_center) / zoom_factor])
        ax.set_ylim([y_center - (y_center - y_lim[0]) / zoom_factor,
                     y_center + (y_lim[1] - y_center) / zoom_factor])


