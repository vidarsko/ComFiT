import numpy as np
import scipy as sp
from typing import Optional
from comfit.tool.tool_complete_field import tool_complete_field
from scipy.interpolate import griddata
import plotly.graph_objects as go
from skimage.measure import marching_cubes
from comfit.tool.tool_set_plot_axis_properties_plotly import tool_set_plot_axis_properties_plotly

def plot_field_in_plane_plotly(
        self,
        field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        **kwargs
    ):
    
    """Plots the field in a plane perpendicular to the given normal vector
    
    Uses scipy.interpolate.griddata and plt.plot_trisurf.

    Args:
        field (array-like): The field to be plotted.
        normal_vector (array-like, optional): The normal vector of the plane. Default is [0,1,0].
        position (array-like, optional): The position of the plane. Default is the middle of the system.
        **kwargs: Keyword arguments for the plot, see https://comfitlib.com/Plotting/. 
    
    Returns:
        The axes containing the plot. (matplotlib.axes.Axes)
    """

    if self.dim != 3:
        raise Exception("The plot in plane function is only defined for 3D fields.")

    kwargs['plot_is_3D'] = True

    # Extend the field if not a complete array is given
    field = tool_complete_field(self, field)

    # Check if the vector field is complex
    if np.iscomplexobj(field):
        print("\033[91mWarning: the provided field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
        print('Max imaginary part: ', np.max(np.imag(field)))
        field = np.real(field)

    fig = kwargs.get('fig', go.Figure())
    colorbar = kwargs.get('colorbar', True)

    # Default values of position and normal vector
    if position is None:
        position = self.rmid

    if normal_vector is None:
        normal_vector = [0, 1, 0]

    normal_vector = np.array(normal_vector) / np.linalg.norm(normal_vector)
    height_above_plane = (self.x - position[0]) * normal_vector[0] + (self.y - position[1]) * normal_vector[1] + (self.z - position[2]) * normal_vector[2]

    verts, faces, _, _ = marching_cubes(height_above_plane, 0)

    # Calculate the centroids of each triangle
    centroids = np.mean(verts[faces], axis=1)

    # Assuming field is defined on the same grid as height_above_plane
    x, y, z = np.mgrid[0:height_above_plane.shape[0], 0:height_above_plane.shape[1], 0:height_above_plane.shape[2]]

    # Flatten the grid for interpolation
    points = np.c_[x.ravel(), y.ravel(), z.ravel()]
    field_values = field.ravel()

    # Interpolate field at the vertices positions
    field_verts = griddata(points, field_values, centroids, method='nearest')

    # Add trace
    fig.add_trace(go.Mesh3d(
        x=(self.xmin + verts[:, 0] * self.dx) / self.a0,
        y=(self.ymin + verts[:, 1] * self.dy) / self.a0,
        z=(self.zmin + verts[:, 2] * self.dz) / self.a0,
        i=faces[:, 0],
        j=faces[:, 1],
        k=faces[:, 2],
        intensity=field_verts,  
        intensitymode='cell',  
        colorscale='Viridis',
        showscale=colorbar
    ))

    kwargs['fig'] = fig
    tool_set_plot_axis_properties_plotly(self, **kwargs)
    return fig
