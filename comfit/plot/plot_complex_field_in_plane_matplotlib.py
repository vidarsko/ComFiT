# Typing imports
from typing import TYPE_CHECKING, Optional, Any, Tuple
if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

# General packages
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from skimage.measure import marching_cubes

# Comfit packages
from comfit.tool import tool_complete_field, \
                        tool_set_plot_axis_properties_matplotlib, \
                        tool_colormap, \
                        tool_matplotlib_define_2D_plot_ax, \
                        tool_matplotlib_define_3D_plot_ax

def plot_complex_field_in_plane_matplotlib(
        self: 'BaseSystem',
        complex_field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        **kwargs: Any
        ) -> Tuple[Figure, Axes]:
    """Plot the complex field in a plane perpendicular to the given normal vector.

    Parameters
    ----------
    self : BaseSystem
        A BaseSystem (or derived) instance.
    complex_field : np.ndarray
        The complex field to be plotted.
    normal_vector : np.ndarray, optional
        The normal vector of the plane. Defaults to [0,1,0].
    position : np.ndarray, optional
        The position of the plane. Defaults to the middle of the system.
    \*\*kwargs : Any
        Keyword arguments for the plot, see https://comfitlib.com/Plotting/.

    Returns
    -------
    tuple[Figure, Axes]
        The figure and axes containing the plot.

    Raises
    ------
    Exception
        If the field dimension is not 3D.
    """


    if self.dim != 3:
        raise Exception("The plot in plane function is only defined for 3D fields.")

    complex_field, fig, ax, kwargs = self.plot_prepare(complex_field, field_type = 'complex', **kwargs)

    ax = tool_matplotlib_define_3D_plot_ax(fig, ax)

    kwargs['plot_is_3D'] = True

    # Default values of position and normal vector
    if position is None:
        position = self.rmid

    if normal_vector is None:
        normal_vector = [0,1,0]

    # Calculate the magnitude and phase of the complex field
    rho = np.abs(complex_field)
    theta = np.angle(complex_field)

    normal_vector = np.array(normal_vector)/np.linalg.norm(normal_vector)
    height_above_plane = (self.x-position[0])*normal_vector[0] + (self.y-position[1])*normal_vector[1] + (self.z-position[2])*normal_vector[2]

    verts, faces, _, _ = marching_cubes(height_above_plane, 0)

    # Calculate the centroids of each triangle
    centroids = np.mean(verts[faces], axis=1)

    # Assuming field is defined on the same grid as height_above_plane
    x, y, z = np.mgrid[0:height_above_plane.shape[0], 0:height_above_plane.shape[1], 0:height_above_plane.shape[2]]

    # Flatten the grid for interpolation
    points = np.c_[x.ravel(), y.ravel(), z.ravel()]
    theta_values = theta.ravel()

    # Interpolate field at the vertices positions
    theta_verts = sp.interpolate.griddata(points, theta_values, centroids, method='nearest')
    rho_verts = sp.interpolate.griddata(points, rho.ravel(), centroids, method='nearest')

    # Normalize field values for color mapping
    theta_normalized = (theta_verts+np.pi) / (2*np.pi)

    # Map normalized field values to colors
    colormap = kwargs['colormap_object']
    colors = colormap(theta_normalized)

    # Blend the colors with white according to rho (normalized)
    colors[:,3] = (rho_verts/np.max(rho_verts)).ravel()

    ax.plot_trisurf((self.xmin+verts[:, 0]*self.dx)/self.a0,
                    (self.ymin+verts[:, 1]*self.dy)/self.a0,
                    faces,
                    (self.zmin+verts[:, 2]*self.dz)/self.a0,
                    facecolor=colors, antialiased=True)


    # Create a colorbar
    if kwargs['colorbar']:
        padding=0.2
        mappable = plt.cm.ScalarMappable(cmap=kwargs['colormap_object'])
        mappable.set_array([])
        mappable.set_clim(-np.pi, np.pi)
        cbar = plt.colorbar(mappable, ax=ax, pad=padding)

        cticks = kwargs.get('cticks', [-np.pi, -2*np.pi/3, -np.pi/3, 0, np.pi/3, 2*np.pi/3, np.pi])
        cbar.set_ticks(cticks)

        cticklabelse = kwargs.get('cticklabels', [r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$'])
        cbar.set_ticklabels(cticklabelse)

    kwargs['grid'] = kwargs.get('grid', True)
    kwargs['ax'] = ax
    tool_set_plot_axis_properties_matplotlib(self, **kwargs)

    return fig, ax