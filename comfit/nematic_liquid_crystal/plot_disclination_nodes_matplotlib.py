import numpy as np
import matplotlib.pyplot as plt
from comfit.tool.tool_complete_field import tool_complete_field
from comfit.tool.tool_set_plot_axis_properties_matplotlib import tool_set_plot_axis_properties_matplotlib

def plot_disclination_nodes_matplotlib(self, disclination_nodes, **kwargs):
    """Plots the discliation nodes in the given axes.

    Args:
        disclination_nodes: A list of dictionaries representing the disclination nodes. Each dictionary contains the following keys:
                                - 'position': The position of the disclination node as a list [x, y].
                                - 'position_index': The index of the position
                                - 'charge' (2D only): The charge of the disclination node.
                                - 'velocity' (currently 2D only): The velocity of the disclination
                                - 'polarization' (2D only): The polarization of the +1/2 disclinations
                                - 'Tangent_vector' (3D only): the tangent of the disclination line
                                - 'Rotation_vector' (3D only): the rotation vector of the disclination line
        -**kwargs: Keyword arguments for the plot.
            See https://comfitlib.com/ClassBaseSystem/
            for a full list of keyword arguments.

    Returns:
        The axes on which the disclination nodes are plotted. (matplotlib.axes.Axes: )
    """

    # Check if an axis object is provided
    fig = kwargs.get('fig', plt.gcf())
    ax = kwargs.get('ax', None)

    if self.dim == 2:

        if ax == None:
            fig.clf()
            ax = plt.gcf().add_subplot(111)

        x_coords_pos = []
        y_coords_pos = []

        x_coords_neg = []
        y_coords_neg = []

        vx_coords_pos = []
        vy_coords_pos = []


        vx_coords_neg = []
        vy_coords_neg = []

        px_coords_pos = []
        py_coords_pos = []



        for disclination in disclination_nodes:

            if disclination['charge'] > 0:
                x_coords_pos.append(disclination['position'][0])
                y_coords_pos.append(disclination['position'][1])
                vx_coords_pos.append(disclination['velocity'][0])
                vy_coords_pos.append(disclination['velocity'][1])
                px_coords_pos.append(disclination['polarization'][0])
                py_coords_pos.append(disclination['polarization'][1])
            else:
                x_coords_neg.append(disclination['position'][0])
                y_coords_neg.append(disclination['position'][1])
                vx_coords_neg.append(disclination['velocity'][0])
                vy_coords_neg.append(disclination['velocity'][1])


        # print(x_coords_pos,y_coords_pos)
        # print(x_coords_neg,y_coords_neg)
        ax.scatter(x_coords_pos, y_coords_pos, marker='+', color='red')
        ax.scatter(x_coords_neg, y_coords_neg, marker='*', color='blue')
        ax.quiver(x_coords_pos, y_coords_pos, vx_coords_pos, vy_coords_pos, color='black')
        ax.quiver(x_coords_neg, y_coords_neg, vx_coords_neg, vy_coords_neg, color='black')
        ax.quiver(x_coords_pos, y_coords_pos, px_coords_pos, py_coords_pos, color='red')
        ax.set_aspect('equal')
        ax.set_facecolor('none')

        ax.set_xlabel('$x/a_0$')
        ax.set_ylabel('$y/a_0$')

        ax.set_xlim([0, self.xmax-self.dx])
        ax.set_ylim([0, self.ymax-self.dy])
        return ax

    elif self.dim == 3:
        # Plotting options
        quiver_scale = 2  # The scale of the quiver arrows

        if ax == None:
            plt.clf()
            ax = plt.gcf().add_subplot(111, projection='3d')
        x_coords = []
        y_coords = []
        z_coords = []

        tx = []
        ty = []
        tz = []

        Ox = []
        Oy = []
        Oz = []

        for disclination in disclination_nodes:
            x_coords.append(disclination['position'][0])
            y_coords.append(disclination['position'][1])
            z_coords.append(disclination['position'][2])

            tx.append(disclination['Tangent_vector'][0])
            ty.append(disclination['Tangent_vector'][1])
            tz.append(disclination['Tangent_vector'][2])

            Ox.append(disclination['Rotation_vector'][0])
            Oy.append(disclination['Rotation_vector'][1])
            Oz.append(disclination['Rotation_vector'][2])

        tx = np.array(tx)
        ty = np.array(ty)
        tz = np.array(tz)

        Ox = np.array(Ox)
        Oy = np.array(Oy)
        Oz = np.array(Oz)


        # ax.scatter(x_coords, y_coords, z_coords, marker='o', color='black')
        ax.quiver(x_coords, y_coords, z_coords, quiver_scale * tx, quiver_scale * ty, quiver_scale * tz,
                    color='blue')
        ax.quiver(x_coords, y_coords, z_coords, quiver_scale * Ox*0.75 , quiver_scale * Oy*0.75 ,
                    quiver_scale * Oz*0.75, color='green')

        ax.set_xlabel('$x/a_0$')
        ax.set_ylabel('$y/a_0$')
        ax.set_zlabel('$z/a_0$')

        ax.set_xlim([0, self.xmax - self.dx])
        ax.set_ylim([0, self.ymax - self.dy])
        ax.set_zlim([0, self.zmax - self.dz])

        ax.set_aspect('equal')
        ax.grid(True)

        return ax