import numpy as np
import matplotlib.pyplot as plt
from comfit.tool.tool_set_plot_axis_properties_matplotlib import tool_set_plot_axis_properties_matplotlib
from comfit.tool import tool_extract_node_arrays


def plot_nodes_matplotlib(self, nodes, **kwargs):
    """
    Plots the nodes.

    Args:
        nodes: The nodes to plot.
        **kwargs: Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/ for a full list of keyword arguments.  
            
    Returns: 
        The axes containing the plot. (matplotlib.axes.Axes)
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
        quiver_scale = 2 # The scale of the quiver arrows

        x_coords = node_arrays['x_coordinates']
        y_coords = node_arrays['y_coordinates']
        z_coords = node_arrays['z_coordinates']

        if node_arrays['tangent_vector_given']:
            tx = np.array(node_arrays['tangent_vector_x_coordinates'])
            ty = np.array(node_arrays['tangent_vector_y_coordinates'])
            tz = np.array(node_arrays['tangent_vector_z_coordinates'])
        
        if node_arrays['velocity_given']:
            vx = np.array(node_arrays['velocity_x_coordinates'])
            vy = np.array(node_arrays['velocity_y_coordinates'])
            vz = np.array(node_arrays['velocity_z_coordinates'])

        if node_arrays['Burgers_vector_given']:
            Bx = np.array(node_arrays['Burgers_vector_x_coordinates'])
            By = np.array(node_arrays['Burgers_vector_y_coordinates'])
            Bz = np.array(node_arrays['Burgers_vector_z_coordinates'])

        if not len(Bx) == 0:
            B2 =Bx**2 + By**2 + Bz**2
            B_norm = np.sqrt(max(B2))
        else:
            B_norm = 1

        ax.scatter(x_coords, y_coords, z_coords, marker='o', color='black')
        if tangent_vector_given:
            ax.quiver(x_coords, y_coords, z_coords, quiver_scale*tx, quiver_scale*ty, quiver_scale*tz, color='blue')
        
        if velocity_given:
            ax.quiver(x_coords, y_coords, z_coords, quiver_scale*vx/v_norm, quiver_scale*vy/v_norm, quiver_scale*vz/v_norm, color='green')
        
        if Burgers_vector_given:
            ax.quiver(x_coords, y_coords, z_coords, quiver_scale*Bx/B_norm, quiver_scale*By/B_norm, quiver_scale*Bz/B_norm, color='red')

    kwargs['fig'] = fig
    kwargs['ax'] = ax
    tool_set_plot_axis_properties_matplotlib(self, **kwargs)

    return fig, ax