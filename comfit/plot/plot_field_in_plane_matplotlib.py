from typing import Optional

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy as sp
from skimage.measure import marching_cubes

from comfit.tool.tool_complete_field import tool_complete_field
from comfit.tool.tool_set_plot_axis_properties_matplotlib import tool_set_plot_axis_properties_matplotlib
from comfit.tool.tool_colormap import tool_colormap
from comfit.tool import tool_matplotlib_define_3D_plot_ax

from mpl_toolkits.mplot3d import Axes3D

def plot_field_in_plane_matplotlib(
        self,
        field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        **kwargs
    ) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    
    """Plots the field in a plane perpendicular to the given normal vector
    
    Uses scipy.interpolate.griddata and plt.plot_trisurf.

    Args:
        field (array-like): The field to be plotted.
        normal_vector (array-like, optional): The normal vector of the plane. Default is [0,1,0].
        position (array-like, optional): The position of the plane. Default is the middle of the system.
        **kwargs: Keyword arguments for the plot, see https://comfitlib.com/Plotting/. 
    
    Returns:
        tuple: A tuple containing:
            - matplotlib.figure.Figure: The figure containing the plot.
            - matplotlib.axes.Axes: The axes containing the plot
    """

    # print('\033[93mWarning: The plot_in_plane function is not yet supported with plotly.\033[0m')

    if self.dim != 3:
        raise Exception("The plot in plane function is only defined for 3D fields.")

    field, fig, ax, kwargs = self.plot_prepare(field, field_type = 'real', **kwargs)

    ax = tool_matplotlib_define_3D_plot_ax(fig, ax) 

    # Default values of position and normal vector
    if position is None:
        position = self.rmid

    if normal_vector is None:
        normal_vector=[0,1,0]

    if ax is None:
        ax = fig.add_subplot(111, projection='3d')

    normal_vector = np.array(normal_vector)/np.linalg.norm(normal_vector)
    height_above_plane = (self.x-position[0])*normal_vector[0] + (self.y-position[1])*normal_vector[1] + (self.z-position[2])*normal_vector[2]

    verts, faces, _, _ = marching_cubes(height_above_plane, 0)

    # Calculate the centroids of each triangle
    centroids = np.mean(verts[faces], axis=1)

    # Assuming field is defined on the same grid as height_above_plane
    x, y, z = np.mgrid[0:height_above_plane.shape[0], 0:height_above_plane.shape[1], 0:height_above_plane.shape[2]]

    # Flatten the grid for interpolation
    points = np.c_[x.ravel(), y.ravel(), z.ravel()]
    field_values = field.ravel()

    # Interpolate field at the vertices positions
    field_verts = sp.interpolate.griddata(points, field_values, centroids, method='nearest')

    # Normalize field values for color mapping
    field_normalized = (field_verts - np.min(field_verts)) / (np.max(field_verts) - np.min(field_verts))

    # Map normalized field values to colors
    colors = kwargs['colormap_object'](field_normalized)

    ax.plot_trisurf((self.xmin+verts[:, 0]*self.dx)/self.a0,
                    (self.ymin+verts[:, 1]*self.dy)/self.a0,
                    faces,
                    (self.zmin+verts[:, 2]*self.dz)/self.a0,
                    facecolor=colors, antialiased=False)

    if kwargs['colorbar']:
            sm = plt.cm.ScalarMappable(cmap=kwargs['colormap_object'])
            sm.set_clim(kwargs['vmin'], kwargs['vmax'])
            plt.colorbar(sm, ax=ax, pad=0.2)

    kwargs['ax'] = ax
    ax = tool_set_plot_axis_properties_matplotlib(self, **kwargs)
    return fig, ax