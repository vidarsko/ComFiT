import numpy as np
import plotly.graph_objects as go
from comfit.tool.tool_complete_field import tool_complete_field
from comfit.tool.tool_set_plot_axis_properties_plotly import tool_set_plot_axis_properties_plotly

def plot_disclination_nodes_plotly(self, disclination_nodes, **kwargs):
    """Plots the disclination nodes in a 2D or 3D Plotly figure.

    Args:
        disclination_nodes: A list of dictionaries representing the disclination nodes. Each dictionary contains:
                                - 'position': The position of the disclination node as a list [x, y] (or [x, y, z] in 3D).
                                - 'charge' (2D only): The charge of the disclination node.
                                - 'velocity' (2D only): The velocity of the disclination node.
                                - 'polarization' (2D only): The polarization of the +1/2 disclinations.
                                - 'Tangent_vector' (3D only): The tangent of the disclination line.
                                - 'Rotation_vector' (3D only): The rotation vector of the disclination line.
        **kwargs: Additional keyword arguments for the plot.

    Returns:
        A Plotly figure (go.Figure).
    """
    fig = kwargs.get('fig', go.Figure())

    if self.dim == 2:
        x_pos, y_pos, vx_pos, vy_pos, px_pos, py_pos = [], [], [], [], [], []
        x_neg, y_neg, vx_neg, vy_neg = [], [], [], []

        for disclination in disclination_nodes:
            if disclination['charge'] > 0:
                x_pos.append(disclination['position'][0])
                y_pos.append(disclination['position'][1])
                vx_pos.append(disclination['velocity'][0])
                vy_pos.append(disclination['velocity'][1])
                px_pos.append(disclination['polarization'][0])
                py_pos.append(disclination['polarization'][1])
            else:
                x_neg.append(disclination['position'][0])
                y_neg.append(disclination['position'][1])
                vx_neg.append(disclination['velocity'][0])
                vy_neg.append(disclination['velocity'][1])

        # Scatter plot for positive and negative charges
        fig.add_trace(go.Scatter(
            x=x_pos, y=y_pos, mode='markers', marker=dict(symbol='cross', color='red', size=10),
            name='+ Charge'
        ))
        fig.add_trace(go.Scatter(
            x=x_neg, y=y_neg, mode='markers', marker=dict(symbol='star', color='blue', size=10),
            name='- Charge'
        ))

        # Quiver plot for velocity and polarization
        fig.add_trace(go.Scatter(
            x=x_pos, y=y_pos, mode='lines+markers',
            marker=dict(size=1),
            line=dict(color='black', width=1),
            name='Velocity (+)',
            customdata=np.column_stack((vx_pos, vy_pos)),
            hovertemplate="Velocity: (%{customdata[0]:.2f}, %{customdata[1]:.2f})"
        ))
        fig.add_trace(go.Scatter(
            x=x_pos, y=y_pos, mode='lines+markers',
            marker=dict(size=1),
            line=dict(color='red', width=1),
            name='Polarization (+)',
            customdata=np.column_stack((px_pos, py_pos)),
            hovertemplate="Polarization: (%{customdata[0]:.2f}, %{customdata[1]:.2f})"
        ))

    elif self.dim == 3:
        x_coords, y_coords, z_coords = [], [], []
        tx, ty, tz = [], [], []
        ox, oy, oz = [], [], []

        for disclination in disclination_nodes:
            x_coords.append(disclination['position'][0])
            y_coords.append(disclination['position'][1])
            z_coords.append(disclination['position'][2])
            tx.append(disclination['Tangent_vector'][0])
            ty.append(disclination['Tangent_vector'][1])
            tz.append(disclination['Tangent_vector'][2])
            ox.append(disclination['Rotation_vector'][0])
            oy.append(disclination['Rotation_vector'][1])
            oz.append(disclination['Rotation_vector'][2])

        # Quiver plot for tangent and rotation vectors
        fig.add_trace(go.Cone(
            x=x_coords, y=y_coords, z=z_coords,
            u=tx, v=ty, w=tz,
            sizemode="scaled", sizeref=0.5,
            anchor="tail",
            colorscale="Blues", name='Tangent Vector'
        ))
        fig.add_trace(go.Cone(
            x=x_coords, y=y_coords, z=z_coords,
            u=ox, v=oy, w=oz,
            sizemode="scaled", sizeref=0.5,
            anchor="tail",
            colorscale="Greens", name='Rotation Vector'
        ))

    # Set axis properties
    kwargs['fig'] = fig
    tool_set_plot_axis_properties_plotly(self, **kwargs)
    return fig
