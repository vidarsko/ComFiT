from typing import TYPE_CHECKING, Any, List, Dict, Tuple

if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

import numpy as np

import plotly.graph_objects as go
import plotly.figure_factory as ff

from comfit.tool import (
    tool_set_plot_axis_properties_plotly,
    tool_plotly_define_2D_plot_ax,
    tool_plotly_define_3D_plot_ax,
    tool_extract_node_arrays
)

def plot_nodes_plotly(
    self: 'BaseSystem',
    nodes: List[Dict],
    **kwargs: Any
) -> Tuple[go.Figure, Dict]:
    """Plot the nodes in the system.

    Parameters
    ----------
    self : BaseSystem
        A BaseSystem (or derived) instance.
    nodes : List[Dict]
        A list of dictionaries representing nodes. Each dictionary must contain:
            - 'position': The position of the node in the defect density array
        Optional keys:
            - 'velocity': The velocity of the node
            - 'charge': The charge of the node
            - 'tangent_vector': The tangent of the node (3D only)
            - 'Burgers_vector': The Burgers vector of the node (for dislocations)
    kwargs : Any
        Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/
        for a full list of keyword arguments.

    Returns
    -------
    Tuple[go.Figure, Dict]
        The figure object and axes dictionary containing the plot configuration.
    """

    # Check if an axis object is provided
    fig = kwargs.get('fig', go.Figure())
    ax = kwargs.get('ax', {'row': 1, 'col': 1, 'nrows': 1, 'ncols': 1})

    # Check if there are nodes to be plotted
    if not nodes:
        return fig, ax

    # Extract node arrays
    node_arrays = tool_extract_node_arrays(self, nodes)

    if self.dim == 2:

        ax = tool_plotly_define_2D_plot_ax(fig, ax)

        if node_arrays['charge_given']:
            if len(node_arrays['x_coordinates_negative']) > 0:
                fig.add_trace(go.Scatter(x=np.array(node_arrays['x_coordinates_negative'])/self.a0, 
                                        y=np.array(node_arrays['y_coordinates_negative'])/self.a0, 
                                        mode='markers', 
                                        marker=dict(symbol='star', color='blue'),
                                        showlegend=False,
                                        xaxis=ax['xN'],
                                        yaxis=ax['yN']))

            if len(node_arrays['x_coordinates_positive']) > 0:
                fig.add_trace(go.Scatter(x=np.array(node_arrays['x_coordinates_positive'])/self.a0, 
                                        y=np.array(node_arrays['y_coordinates_positive'])/self.a0, 
                                        mode='markers', 
                                        marker=dict(symbol='cross', color='red'),
                                        showlegend=False, 
                                        xaxis=ax['xN'],
                                        yaxis=ax['yN']))
        
        else: # If no charge is provided, plot all nodes as positive
            fig.add_trace(go.Scatter(x=np.array(node_arrays['x_coordinates'])/self.a0, 
                                    y=np.array(node_arrays['y_coordinates'])/self.a0, 
                                    mode='markers', 
                                    marker=dict(symbol='circle', color='black'),
                                    showlegend=False,
                                    xaxis=ax['xN'],
                                    yaxis=ax['yN']))
        
        # Velocities
        if node_arrays['velocity_given']:
            
            vx = np.array(node_arrays['velocity_x_coordinates'])
            vy = np.array(node_arrays['velocity_y_coordinates'])

            velocity_scale = 0.1*min(self.size_x, self.size_y)/self.a0

            vnorm = np.sqrt(np.max(vx**2 + vy**2))
                
            vx = velocity_scale*vx/vnorm
            vy = velocity_scale*vy/vnorm

            fig1 = ff.create_quiver(x=np.array(node_arrays['x_coordinates'])/self.a0, 
                                    y=np.array(node_arrays['y_coordinates'])/self.a0, 
                                    u=vx, 
                                    v=vy, 
                                    line=dict(width=1, color='black'),
                                    scale=1,
                                    showlegend=False,
                                    xaxis=ax['xN'],
                                    yaxis=ax['yN'])
                
            fig.add_traces(data=fig1.data)

        # Burgers vectors
        if node_arrays['Burgers_vector_given']:

            bvx = np.array(node_arrays['Burgers_vector_x_coordinates'])
            bvy = np.array(node_arrays['Burgers_vector_y_coordinates'])

            Burgers_vector_scale = 0.1*min(self.size_x, self.size_y)/self.a0

            bnorm = np.sqrt(np.max(bvx**2 + bvy**2))
                
            bvx = Burgers_vector_scale*bvx/bnorm
            bvy = Burgers_vector_scale*bvy/bnorm

            fig1 = ff.create_quiver(x=np.array(node_arrays['x_coordinates'])/self.a0, 
                                    y=np.array(node_arrays['y_coordinates'])/self.a0, 
                                    u=bvx, 
                                    v=bvy, 
                                    line=dict(width=1, color='red'),
                                    scale=1,
                                    showlegend=False,
                                    xaxis=ax['xN'],
                                    yaxis=ax['yN'])
                
            fig.add_traces(data=fig1.data)

    elif self.dim == 3:

        ax = tool_plotly_define_3D_plot_ax(fig, ax)

        x = np.array(node_arrays['x_coordinates'])
        y = np.array(node_arrays['y_coordinates'])
        z = np.array(node_arrays['z_coordinates'])

        # Plotting options
        quiver_scale = 2 # The scale of the quiver arrows
        if node_arrays['tangent_vector_given']:

            tx = np.array(node_arrays['tangent_vector_x_coordinates'])
            ty = np.array(node_arrays['tangent_vector_y_coordinates'])
            tz = np.array(node_arrays['tangent_vector_z_coordinates'])

            fig.add_trace(go.Cone(x=x / self.a0,
                                  y=y / self.a0,
                                  z=z / self.a0,
                                  u=tx,
                                  v=ty,
                                  w=tz,
                                  scene=ax['sceneN'],
                                  colorscale='Blues',
                                  sizemode='scaled',
                                  sizeref=0))

        if node_arrays['velocity_given']:

            vx = np.array(node_arrays['velocity_x_coordinates'])
            vy = np.array(node_arrays['velocity_y_coordinates'])
            vz = np.array(node_arrays['velocity_z_coordinates'])

            v2 =vx**2 + vy**2 + vz**2
            v_norm = np.sqrt(max(v2))

            if v_norm > 0:
                vx /= v_norm
                vy /= v_norm
                vz /= v_norm

            fig.add_trace(go.Cone(x=x / self.a0,
                                  y=y / self.a0,
                                  z=z / self.a0,
                                  u=vx,
                                  v=vy,
                                  w=vz,
                                  scene=ax['sceneN'],
                                  colorscale='Greens',
                                  sizemode='scaled',
                                  sizeref=0))

        if node_arrays['rotation_vector_given']:
            rx = np.array(node_arrays['rotation_vector_x_coordinates'])
            ry = np.array(node_arrays['rotation_vector_y_coordinates'])
            rz = np.array(node_arrays['rotation_vector_z_coordinates'])



            fig.add_trace(go.Cone(x=x / self.a0,
                                  y=y / self.a0,
                                  z=z / self.a0,
                                  u=rx,
                                  v=ry,
                                  w=rz,
                                  scene=ax['sceneN'],
                                  colorscale='Reds',
                                  sizemode='scaled',
                                  sizeref=0))

            #ax.scatter(node_arrays['x_coordinates'], node_arrays['y_coordinates'], z_coords, marker='o', color='black')\
            # fig.add_trace(go.Scatter(x=node_arrays['x_coordinates']/self.a0, 
            #                             y=node_arrays['y_coordinates']/self.a0, 
            #                             mode='markers', 
            #                             marker=dict(symbol='circle', color='black')))   #TODO: Check if this works




        

    kwargs['fig'] = fig
    kwargs['ax'] = ax
    tool_set_plot_axis_properties_plotly(self, **kwargs)

    return fig, ax