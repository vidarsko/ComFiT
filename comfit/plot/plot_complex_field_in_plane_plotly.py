# Typing imports
from typing import TYPE_CHECKING, Optional, Any
if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

# General packages
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from skimage.measure import marching_cubes

# Local application imports
from comfit.tool import tool_complete_field,\
                        tool_colormap, \
                        tool_set_plot_axis_properties_plotly, \
                        tool_plotly_define_3D_plot_ax, \
                        tool_plotly_colorbar

def plot_complex_field_in_plane_plotly(
        self: 'BaseSystem',
        complex_field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        **kwargs: Any
        ) -> tuple[go.Figure, dict]:
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
    tuple[go.Figure, dict]
        The figure and axes dictionary containing the plot of the complex field.

    Raises
    ------
    Exception
        If the dimension of the system is not 3.
    """

    if self.dim != 3:
        raise Exception("The plot in plane function is only defined for 3D fields.")

    kwargs['colormap'] = kwargs.get('colormap', 'angle') # Override the default colormap with 'angle'
    complex_field, fig, ax, kwargs = self.plot_prepare(complex_field, field_type = 'complex', **kwargs)

    ax = tool_plotly_define_3D_plot_ax(fig, ax) #Defines sceneN, plot_dimension

    # Default values of position and normal vector
    if position is None:
        position = self.rmid

    if normal_vector is None:
        normal_vector = [0,1,0]

    # Calculate the magnitude and phase of the complex field
    rho = np.abs(complex_field)
    theta = np.angle(complex_field)

    reals = np.real(complex_field)
    imags = np.imag(complex_field)

    # Kewyord arguments
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
    rho_values = rho.ravel()

    # Interpolate field at the vertices positions
    # theta_verts = sp.interpolate.griddata(points, theta_values, centroids, method='nearest')
    # rho_verts = sp.interpolate.griddata(points, rho_values, centroids, method='nearest')

    interpolation_method = kwargs.get('interpolation_method', 'nearest') #TODO: This line is inconsistent in naming with
    print("Interpolating points with method: ' ", interpolation_method, "'.")
    print("If this process is slow, consider passing 'interpolation_method='nearest' with the plot_complex_field_in_plane function.") #TODO: this line. (interpolation vs interpolation_method)
    print("That will speed up the process, but the plot may look less smooth.")
    reals_verts = sp.interpolate.griddata(points, reals.ravel(), centroids, method=interpolation_method)
    imags_verts = sp.interpolate.griddata(points, imags.ravel(), centroids, method=interpolation_method)
    print("Interpolation done.")

    theta_verts = np.arctan2(imags_verts, reals_verts)
    rho_verts = np.abs(reals_verts + 1j*imags_verts)

    # Normalize field values for color mapping
    theta_normalized = (theta_verts+np.pi) / (2*np.pi)

    # Map normalized field values to colors
    plt_colormap_object = tool_colormap(kwargs['colormap'], plot_lib='matplotlib')
    colors = plt_colormap_object(theta_normalized)

    # Blend the colors with white according to rho (normalized)
    colors[:,3] = (rho_verts/np.max(rho_verts)).ravel()

    # Create the mesh object
    x_new = (verts[:, 0] * self.dx + self.xmin) / self.a0
    y_new = (verts[:, 1] * self.dy + self.ymin) / self.a0
    z_new = (verts[:, 2] * self.dz + self.zmin) / self.a0

    mesh = go.Mesh3d(
                x=x_new, 
                y=y_new, 
                z=z_new,
                i=faces[:, 0], 
                j=faces[:, 1], 
                k=faces[:, 2],
                customdata=np.array([theta_verts/np.pi, rho_verts]).T,
                facecolor=colors,  # Set color for each face
                showscale=False,
                hovertemplate=kwargs['xlabel']+': %{x:.2f}<br>'+\
                                kwargs['ylabel']+': %{y:.2f}<br>'+\
                                kwargs['zlabel']+': %{z:.2f}<br>'+\
                                'amplitude: %{customdata[1]:.2e}<br>'+\
                                'phase: %{customdata[0]:.2f} π',
                name = '',
                scene=ax['sceneN']
            )
    fig.add_trace(mesh)

    # Create a colorbar
    # if colorbar:
    #     mappable = plt.cm.ScalarMappable(cmap=tool_colormap_angle())
    #     mappable.set_array([])
    #     mappable.set_clim(-np.pi, np.pi)
    #     cbar = plt.colorbar(mappable, ax=ax, pad=0.2)
    #     cbar.set_ticks(np.array([-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi]))
    #     cbar.set_ticklabels([r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$'])


    if kwargs['colorbar'] and not(ax['colorbar']) and not(kwargs['field_is_nan']):
        ax['colormap_object'] = kwargs['colormap_object']
        fig.add_trace(tool_plotly_colorbar(ax, type='angle'))
        ax['colorbar'] = True

    kwargs['fig'] = fig
    kwargs['ax'] = ax
    tool_set_plot_axis_properties_plotly(self, **kwargs)

    return fig, ax