import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def plot_vortex_nodes_matplotlib(self, vortex_nodes: list[dict], **kwargs):
        """Plots the vortex nodes in the system.

        Args:
            vortex_nodes (list): A list of dictionaries representing the vortex nodes. Each dictionary contains the following keys:
                                 - 'position_index': The position index of the vortex node in the defect density array.
                                 - 'charge': The charge of the vortex node.
                                 - 'position': The position of the vortex node as a list [x, y].
                                 - 'velocity': The velocity of the vortex node as a list [vx, vy].
            -**kwargs: Keyword arguments for the plot.
                See https://comfitlib.com/ClassBaseSystem/
                for a full list of keyword arguments.
        
        Returns:
            The axes on which the vortex nodes are plotted. (matplotlib.axes.Axes: )
        """
        # Check if an axis object is provided
        fig = kwargs.get('fig', plt.gcf())
        ax = kwargs.get('ax', None)
        if self.dim == 2:

            if ax == None:
                ax = plt.gcf().add_subplot(111)

            x_coords_pos = []
            y_coords_pos = []

            x_coords_neg = []
            y_coords_neg = []

            vx_coords_pos = []
            vy_coords_pos = []

            vx_coords_neg = []
            vy_coords_neg = []

            for vortex in vortex_nodes:

                if vortex['charge'] > 0:
                    x_coords_pos.append(vortex['position'][0])
                    y_coords_pos.append(vortex['position'][1])
                    vx_coords_pos.append(vortex['velocity'][0])
                    vy_coords_pos.append(vortex['velocity'][1])
                else:
                    x_coords_neg.append(vortex['position'][0])
                    y_coords_neg.append(vortex['position'][1])
                    vx_coords_neg.append(vortex['velocity'][0])
                    vy_coords_neg.append(vortex['velocity'][1])

            # print(x_coords_pos,y_coords_pos)
            # print(x_coords_neg,y_coords_neg)
            ax.scatter(x_coords_pos, y_coords_pos, marker='+', color='red')
            ax.scatter(x_coords_neg, y_coords_neg, marker='*', color='blue')
            ax.quiver(x_coords_pos, y_coords_pos, vx_coords_pos, vy_coords_pos, color='black')
            ax.quiver(x_coords_neg, y_coords_neg, vx_coords_neg, vy_coords_neg, color='black')
            ax.set_aspect('equal')
            ax.set_facecolor('none')

            ax.set_xlabel('$x/a_0$')
            ax.set_ylabel('$y/a_0$')

            ax.set_xlim([0, self.xmax-self.dx])
            ax.set_ylim([0, self.ymax-self.dy])

        elif self.dim == 3:
            # Plotting options
            quiver_scale = 2 # The scale of the quiver arrows

            if ax == None:
                ax = plt.gcf().add_subplot(111, projection='3d')
               # ax = plt.gca()

            x_coords = []
            y_coords = []
            z_coords = []

            tx = []
            ty = []
            tz = []

            vx = []
            vy = []
            vz = []


            for vortex in vortex_nodes:
                x_coords.append(vortex['position'][0])
                y_coords.append(vortex['position'][1])
                z_coords.append(vortex['position'][2])

                tx.append(vortex['tangent_vector'][0])
                ty.append(vortex['tangent_vector'][1])
                tz.append(vortex['tangent_vector'][2])

                vx.append(vortex['velocity'][0])
                vy.append(vortex['velocity'][1])
                vz.append(vortex['velocity'][2])

            tx = np.array(tx)
            ty = np.array(ty)
            tz = np.array(tz)

            vx = np.array(vx)
            vy = np.array(vy)
            vz = np.array(vz)

            if not len(vx) == 0:
                v2 =vx**2 + vy**2 + vz**2
                v_norm = np.sqrt(max(v2))
            else:
                v_norm = 1

            #ax.scatter(x_coords, y_coords, z_coords, marker='o', color='black')
            ax.quiver(x_coords, y_coords, z_coords, quiver_scale*tx, quiver_scale*ty, quiver_scale*tz, color='blue')
            ax.quiver(x_coords, y_coords, z_coords, quiver_scale*vx/v_norm, quiver_scale*vy/v_norm, quiver_scale*vz/v_norm, color='green')

            ax.set_xlabel('$x/a_0$')
            ax.set_ylabel('$y/a_0$')
            ax.set_zlabel('$z/a_0$')

            ax.set_xlim([0, self.xmax-self.dx])
            ax.set_ylim([0, self.ymax-self.dy])
            ax.set_zlim([0, self.zmax-self.dz])

            ax.set_aspect('equal')
        ax.grid(True)

        return ax