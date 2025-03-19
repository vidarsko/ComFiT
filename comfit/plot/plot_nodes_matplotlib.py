from typing import TYPE_CHECKING, Any, Dict, List, Tuple

if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

import numpy as np
import matplotlib.pyplot as plt

from comfit.tool import (
    tool_set_plot_axis_properties_matplotlib,
    tool_extract_node_arrays,
tool_matplotlib_define_3D_plot_ax
)

def plot_nodes_matplotlib(
        self: 'BaseSystem',
        nodes: List[Dict],
        **kwargs: Any
        ) -> Tuple[plt.Figure, plt.Axes]:
    """Plot the nodes of the system.

    Create a matplotlib plot of the nodes in the system, including their positions,
    charges, velocities and Burgers vectors if available.

    Parameters
    ----------
    self : BaseSystem
        A BaseSystem (or derived) instance.
    nodes : List[Dict]
        List of node dictionaries containing position and property information
    kwargs : Any
        Keyword arguments for customizing the plot. See 
        https://comfitlib.com/ClassBaseSystem/ for full list.

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        The figure and axes objects containing the plot
    """

    # Check if an axis object is provided
    fig = kwargs.get('fig', plt.gcf())
    ax = kwargs.get('ax', None)

    # Check if there are nodes to be plotted, if not return the axes
    if not nodes:
        return fig, ax

    node_arrays = tool_extract_node_arrays(self, nodes)
    
    if self.dim == 2:

        if ax == None:
            fig.clf()
            ax = fig.add_subplot(111)

        x_coords = np.array(node_arrays['x_coordinates'])
        y_coords = np.array(node_arrays['y_coordinates'])

        if node_arrays['charge_given']:
            x_coords_positive = np.array(node_arrays['x_coordinates_positive'])
            y_coords_positive = np.array(node_arrays['y_coordinates_positive'])

            ax.scatter(x_coords_positive/self.a0, y_coords_positive/self.a0, marker='+', color='red')

            x_coords_negative = np.array(node_arrays['x_coordinates_negative'])
            y_coords_negative = np.array(node_arrays['y_coordinates_negative'])

            ax.scatter(x_coords_negative/self.a0, y_coords_negative/self.a0, marker='o', color='blue')
            
        else:
            ax.scatter(x_coords/self.a0, y_coords/self.a0, marker='o', color='black')

        if node_arrays['velocity_given']:
            vx_coords = np.array(node_arrays['velocity_x_coordinates'])
            vy_coords = np.array(node_arrays['velocity_y_coordinates'])
            ax.quiver(x_coords/self.a0, y_coords/self.a0, vx_coords, vy_coords, color='black')

        if node_arrays['Burgers_vector_given']:
            Bx_coords = np.array(node_arrays['Burgers_vector_x_coordinates'])
            By_coords = np.array(node_arrays['Burgers_vector_y_coordinates'])
            ax.quiver(x_coords/self.a0, y_coords/self.a0, Bx_coords, By_coords, color='red')

    elif self.dim == 3:
        # Plotting options

        if ax == None:
            fig.clf()
            ax= tool_matplotlib_define_3D_plot_ax(fig, ax)

        quiver_scale = 2 # The scale of the quiver arrows

        x_coords = node_arrays['x_coordinates']
        y_coords = node_arrays['y_coordinates']
        z_coords = node_arrays['z_coordinates']

        if node_arrays['tangent_vector_given']:
            tx = np.array(node_arrays['tangent_vector_x_coordinates'])
            ty = np.array(node_arrays['tangent_vector_y_coordinates'])
            tz = np.array(node_arrays['tangent_vector_z_coordinates'])

            ax.quiver(x_coords, y_coords, z_coords, quiver_scale * tx, quiver_scale * ty, quiver_scale * tz,
                      color='blue')
        
        if node_arrays['velocity_given']:
            vx = np.array(node_arrays['velocity_x_coordinates'])
            vy = np.array(node_arrays['velocity_y_coordinates'])
            vz = np.array(node_arrays['velocity_z_coordinates'])

            v_norm = np.sqrt(vx**2 +vy**2+vz**2)

            ax.quiver(x_coords, y_coords, z_coords, quiver_scale * vx / v_norm, quiver_scale * vy / v_norm,
                      quiver_scale * vz / v_norm, color='green')

        if node_arrays['Burgers_vector_given']:
            Bx = np.array(node_arrays['Burgers_vector_x_coordinates'])
            By = np.array(node_arrays['Burgers_vector_y_coordinates'])
            Bz = np.array(node_arrays['Burgers_vector_z_coordinates'])

            if not len(Bx) == 0:
                B2 = Bx ** 2 + By ** 2 + Bz ** 2
                B_norm = np.sqrt(max(B2))
            else:
                B_norm = 1

            ax.quiver(x_coords, y_coords, z_coords, quiver_scale * Bx / B_norm, quiver_scale * By / B_norm,
                      quiver_scale * Bz / B_norm, color='red')

        if node_arrays['rotation_vector_given']:
            rx = np.array(node_arrays['rotation_vector_x_coordinates'])
            ry = np.array(node_arrays['rotation_vector_y_coordinates'])
            rz = np.array(node_arrays['rotation_vector_z_coordinates'])

            ax.quiver(x_coords, y_coords, z_coords, quiver_scale * rx, quiver_scale * ry, quiver_scale * rz,
                      color='black')



        ax.scatter(x_coords, y_coords, z_coords, marker='o', color='black')


        



    kwargs['fig'] = fig
    kwargs['ax'] = ax
    tool_set_plot_axis_properties_matplotlib(self, **kwargs)

    return fig, ax