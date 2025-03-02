import numpy as np
import plotly.graph_objects as go
import plotly.figure_factory as ff
from comfit.tool import tool_set_plot_axis_properties_plotly, tool_plotly_define_2D_plot_ax, tool_plotly_define_3D_plot_ax

def plot_nodes_plotly(self, nodes, **kwargs):
        """Plots the  nodes in the system.

        Args:
            nodes (list): A list of dictionaries representing nodes. The dicitonary must contain:
                                 - 'position': The position of the  node in the defect density array.
                                 The following keys are optional:
                                 - 'velocity': The velocity of the node.
                                 - 'charge': The charge of the node.
                                 - 'tangent_vector': The tangent of the node (3D only).
                                 - 'Burgers_vector': The Burgers vector of the node (for dislocations).
            -**kwargs: Keyword arguments for the plot.
                See https://comfitlib.com/ClassBaseSystem/
                for a full list of keyword arguments.
        
        Returns:
            The axes on which the  nodes are plotted. (matplotlib.axes.Axes: )
        """

        # Check if an axis object is provided
        fig = kwargs.get('fig', go.Figure())
        ax = kwargs.get('ax', {'row': 1, 'col': 1, 'nrows': 1, 'ncols': 1})

        # Check if there are nodes to be plotted
        if not nodes:
            return fig, ax

        if self.dim == 2:

            ax = tool_plotly_define_2D_plot_ax(ax, fig)

            # Positions
            x_coords = []
            y_coords = []

            x_coords_pos = []
            y_coords_pos = []

            x_coords_neg = []
            y_coords_neg = []

            charge_given = 'charge' in nodes[0].keys()

            if charge_given:

                for node in nodes:
                    x_coords.append(node['position'][0])
                    y_coords.append(node['position'][1])

                    if node['charge'] > 0:
                        x_coords_pos.append(node['position'][0])
                        y_coords_pos.append(node['position'][1])

                    else:
                        x_coords_neg.append(node['position'][0])
                        y_coords_neg.append(node['position'][1])

            else: # No charge given
                for node in nodes:
                    x_coords_pos.append(node['position'][0])
                    y_coords_pos.append(node['position'][1])
                
                x_coords = x_coords_pos
                y_coords = y_coords_pos

            if charge_given:
                if len(x_coords_neg) > 0:
                    fig.add_trace(go.Scatter(x=np.array(x_coords_neg)/self.a0, 
                                            y=np.array(y_coords_neg)/self.a0, 
                                            mode='markers', 
                                            marker=dict(symbol='star', color='blue'),
                                            showlegend=False,
                                            xaxis=ax['xN'],
                                            yaxis=ax['yN']))

                if len(x_coords_pos) > 0:
                    fig.add_trace(go.Scatter(x=np.array(x_coords_pos)/self.a0, 
                                            y=np.array(y_coords_pos)/self.a0, 
                                            mode='markers', 
                                            marker=dict(symbol='cross', color='red'),
                                            showlegend=False, 
                                            xaxis=ax['xN'],
                                            yaxis=ax['yN']))
            
            else: # If no charge is provided, plot all nodes as positive
                fig.add_trace(go.Scatter(x=np.array(x_coords)/self.a0, 
                                        y=np.array(y_coords)/self.a0, 
                                        mode='markers', 
                                        marker=dict(symbol='circle', color='black'),
                                        showlegend=False,
                                        xaxis=ax['xN'],
                                        yaxis=ax['yN']))
            
            # Velocities
            if 'velocity' in nodes[0].keys():
                velocity_scale = 0.1*min(self.size_x, self.size_y)/self.a0

                vx = []
                vy = []
                
                for node in nodes:
                    vx.append(node['velocity'][0])
                    vy.append(node['velocity'][1])
                
                vx = np.array(vx)
                vy = np.array(vy)

                vnorm = np.sqrt(np.max(vx**2 + vy**2))
                    
                vx = velocity_scale*vx/vnorm
                vy = velocity_scale*vy/vnorm

                fig1 = ff.create_quiver(x=np.array(x_coords)/self.a0, 
                                        y=np.array(y_coords)/self.a0, 
                                        u=vx, 
                                        v=vy, 
                                        line=dict(width=1, color='black'),
                                        scale=1,
                                        showlegend=False,
                                        xaxis=ax['xN'],
                                        yaxis=ax['yN'])
                    
                fig.add_traces(data=fig1.data)

            # Burgers vectors
            if 'Burgers_vector' in nodes[0].keys():
                Burgers_vector_scale = 0.1*min(self.size_x, self.size_y)/self.a0

                bvx = []
                bvy = []
                
                for node in nodes:
                    bvx.append(node['Burgers_vector'][0])
                    bvy.append(node['Burgers_vector'][1])
                
                bvx = np.array(bvx)
                bvy = np.array(bvy)

                bnorm = np.sqrt(np.max(bvx**2 + bvy**2))
                    
                bvx = Burgers_vector_scale*bvx/bnorm
                bvy = Burgers_vector_scale*bvy/bnorm

                fig1 = ff.create_quiver(x=np.array(x_coords)/self.a0, 
                                        y=np.array(y_coords)/self.a0, 
                                        u=bvx, 
                                        v=bvy, 
                                        line=dict(width=1, color='red'),
                                        scale=1,
                                        showlegend=False,
                                        xaxis=ax['xN'],
                                        yaxis=ax['yN'])
                    
                fig.add_traces(data=fig1.data)

        elif self.dim == 3:

            ax = tool_plotly_define_3D_plot_ax(ax, fig)

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


            for node in nodes:
                x_coords.append(node['position'][0])
                y_coords.append(node['position'][1])
                z_coords.append(node['position'][2])

                tx.append(node['tangent_vector'][0])
                ty.append(node['tangent_vector'][1])
                tz.append(node['tangent_vector'][2])

                vx.append(node['velocity'][0])
                vy.append(node['velocity'][1])
                vz.append(node['velocity'][2])

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
                # fig.add_trace(go.Scatter(x=x_coords/self.a0, 
                #                             y=y_coords/self.a0, 
                #                             mode='markers', 
                #                             marker=dict(symbol='circle', color='black')))   #TODO: Check if this works

                fig.add_trace(go.Cone(x=x_coords/self.a0, 
                                        y=y_coords/self.a0, 
                                        z=z_coords/self.a0, 
                                        u=tx, 
                                        v=ty, 
                                        w=tz, 
                                        scene=ax['sceneN'],
                                        colorscale='Blues', 
                                        sizemode='scaled', 
                                        sizeref=0))

                fig.add_trace(go.Cone(x=x_coords/self.a0, 
                                        y=y_coords/self.a0, 
                                        z=z_coords/self.a0, 
                                        u=vx, 
                                        v=vy, 
                                        w=vz, 
                                        scene=ax['sceneN'],
                                        colorscale='Greens', 
                                        sizemode='scaled', 
                                        sizeref=0))
            

        kwargs['fig'] = fig
        tool_set_plot_axis_properties_plotly(self, **kwargs)

        return fig