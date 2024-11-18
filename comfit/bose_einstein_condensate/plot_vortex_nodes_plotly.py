import numpy as np
import plotly.graph_objects as go
import plotly.figure_factory as ff
from comfit.tool import tool_set_plot_axis_properties_plotly

def plot_vortex_nodes_plotly(self, vortex_nodes, **kwargs):
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
        fig = kwargs.get('fig', go.Figure())

        if self.dim == 2:

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

            velociy_scale = 0.1*min(self.size_x, self.size_y)

            vx_coords_neg = np.array(vx_coords_neg)
            vy_coords_neg = np.array(vy_coords_neg)

            if len(vx_coords_neg) > 0:
                vnorm_neg = np.max(np.sqrt(vx_coords_neg**2 + vy_coords_neg**2))
            else: 
                vnorm_neg = 0
    
            vx_coords_pos = np.array(vx_coords_pos)
            vy_coords_pos = np.array(vy_coords_pos)

            if len(vx_coords_pos) > 0:
                vnorm_pos = np.max(np.sqrt(vx_coords_pos**2 + vy_coords_pos**2))
            else:
                vnorm_pos = 0

            vnorm = max(vnorm_neg, vnorm_pos)

            if len(x_coords_neg) > 0:
                fig.add_trace(go.Scatter(x=x_coords_neg, 
                                        y=y_coords_neg, 
                                        mode='markers', 
                                        marker=dict(symbol='star', color='blue'),
                                        showlegend=False))
                
                vx_coords_neg = velociy_scale*vx_coords_neg/vnorm
                vy_coords_neg = velociy_scale*vy_coords_neg/vnorm

                fig1 = ff.create_quiver(x=x_coords_neg, 
                                        y=y_coords_neg, 
                                        u=vx_coords_neg, 
                                        v=vy_coords_neg, 
                                        line=dict(width=1, color='black'), 
                                        scale=1,
                                        showlegend=False)

                fig.add_traces(data=fig1.data)

            if len(x_coords_pos) > 0:
                fig.add_trace(go.Scatter(x=x_coords_pos, 
                                        y=y_coords_pos, 
                                        mode='markers', 
                                        marker=dict(symbol='cross', color='red'),
                                        showlegend=False))

                
                vx_coords_pos = velociy_scale*vx_coords_pos/vnorm
                vy_coords_pos = velociy_scale*vy_coords_pos/vnorm

                fig2 = ff.create_quiver(x=x_coords_pos, 
                                        y=y_coords_pos, 
                                        u=vx_coords_pos, 
                                        v=vy_coords_pos, 
                                        line=dict(width=1, color='black'),
                                        scale=1,
                                        showlegend=False)
                
                fig.add_traces(data=fig2.data)

        elif self.dim == 3:
            # Plotting options
            quiver_scale = 2 # The scale of the quiver arrows

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

            if len(x_coords) > 0:
                #ax.scatter(x_coords, y_coords, z_coords, marker='o', color='black')\
                fig.add_trace(go.Scatter(x=x_coords, y=y_coords, mode='markers', marker=dict(symbol='circle', color='black')))

                fig.add_trace(go.Cone(x=x_coords, y=y_coords, z=z_coords, u=tx, v=ty, w=tz, colorscale='Blues', sizemode='scaled', sizeref=0))

                fig.add_trace(go.Cone(x=x_coords, y=y_coords, z=z_coords, u=vx, v=vy, w=vz, colorscale='Greens', sizemode='scaled', sizeref=0))
            

        kwargs['fig'] = fig
        tool_set_plot_axis_properties_plotly(self, **kwargs)

        return fig