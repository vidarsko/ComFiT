import numpy as np
import matplotlib.pyplot as plt
from comfit.tool.tool_set_plot_axis_properties_matplotlib import tool_set_plot_axis_properties_matplotlib

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

    velocity_given = ('velocity' in nodes[0].keys())
    Burgers_vector_given = ('Burgers_vector' in nodes[0].keys())
    tangent_vector_given = ('tangent_vector' in nodes[0].keys())

    if self.dim == 2:

        if ax == None:
            fig.clf()
            ax = fig.add_subplot(111)

        x_coords = []
        y_coords = []

        if velocity_given:
            vx_coords = []
            vy_coords = []

        if Burgers_vector_given:
            Bx_coords = []
            By_coords = []

        
        for node in nodes:

            x_coords.append(node['position'][0])
            y_coords.append(node['position'][1])

            if velocity_given:
                vx_coords = node['velocity'][0]
                vy_coords = node['velocity'][1]

            if Burgers_vector_given:
                Bx_coords.append(node['Burgers_vector'][0])
                By_coords.append(node['Burgers_vector'][1])

        x_coords = np.array(x_coords)
        y_coords = np.array(y_coords)

        # print(x_coords_pos,y_coords_pos)
        # print(x_coords_neg,y_coords_neg)
        ax.scatter(x_coords/self.a0, y_coords/self.a0, marker='o', color='black')

        if velocity_given:
            ax.quiver(x_coords, y_coords, vx_coords, vy_coords, color='black')

        if Burgers_vector_given:
            ax.quiver(x_coords/self.a0, y_coords/self.a0, Bx_coords, By_coords, color='red')

    elif self.dim == 3:
        # Plotting options
        quiver_scale = 2 # The scale of the quiver arrows

        if ax == None:
            plt.clf()
            ax = plt.gcf().add_subplot(111, projection='3d')

        x_coords = []
        y_coords = []
        z_coords = []

        tx = []
        ty = []
        tz = []

        # vx = []
        # vy = []
        # vz = []

        Bx = []
        By = []
        Bz = []


        for node in nodes:
            x_coords.append(node['position'][0])
            y_coords.append(node['position'][1])
            z_coords.append(node['position'][2])

            if tangent_vector_given:
                tx.append(node['tangent_vector'][0])
                ty.append(node['tangent_vector'][1])
                tz.append(node['tangent_vector'][2])

            if velocity_given:
                vx.append(node['velocity'][0])
                vy.append(node['velocity'][1])
                vz.append(node['velocity'][2])

            if Burgers_vector_given:
                Bx.append(node['Burgers_vector'][0])
                By.append(node['Burgers_vector'][1])
                Bz.append(node['Burgers_vector'][2])

        if tangent_vector_given:
            tx = np.array(tx)
            ty = np.array(ty)
            tz = np.array(tz)

        if velocity_given:
            vx = np.array(vx)
            vy = np.array(vy)
            vz = np.array(vz)

        if Burgers_vector_given:
            Bx = np.array(Bx)
            By = np.array(By)
            Bz = np.array(Bz)

        # if not len(vx) == 0:
        #     v2 =vx**2 + vy**2 + vz**2
        #     v_norm = np.sqrt(max(v2))
        # else:
        #     v_norm = 1

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