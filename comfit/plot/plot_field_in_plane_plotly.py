from typing import TYPE_CHECKING, Optional, Any, Tuple, Dict
if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

# Standard library imports
import numpy as np
from scipy.interpolate import griddata
from skimage.measure import marching_cubes

# Third-party library imports
import plotly.graph_objects as go

# Local imports
from comfit.tool.tool_complete_field import tool_complete_field
from comfit.tool.tool_set_plot_axis_properties_plotly import tool_set_plot_axis_properties_plotly
from comfit.tool import tool_plotly_define_3D_plot_ax, tool_plotly_colorbar

def plot_field_in_plane_plotly(
        self: 'BaseSystem',
        field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        **kwargs: Any
        ) -> Tuple[go.Figure, Dict]:
    """Plot the field in a plane perpendicular to the given normal vector.
    
    Uses scipy.interpolate.griddata and plotly's Mesh3d for visualization.

    Parameters
    ----------
    self : BaseSystem
        A BaseSystem (or derived) instance.
    field : np.ndarray
        The field to be plotted
    normal_vector : np.ndarray, optional
        The normal vector of the plane. Defaults to [0,1,0]
    position : np.ndarray, optional
        The position of the plane. Defaults to the middle of the system
    \*\*kwargs : Any
        Additional keyword arguments for plot customization. See https://comfitlib.com/Plotting/

    Returns
    -------
    Tuple[go.Figure, dict]
        The figure object and axes dictionary containing the plot

    Raises
    ------
    Exception
        If the field dimension is not 3D
    """

    if self.dim != 3:
        raise Exception("The plot in plane function is only defined for 3D fields.")

    field, fig, ax, kwargs = self.plot_prepare(field, field_type = 'real', **kwargs)

    ax = tool_plotly_define_3D_plot_ax(fig, ax) #Defines sceneN, plot_dimension


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

    ax['vmin'] = kwargs.get('vmin', np.min(field_verts))
    ax['vmax'] = kwargs.get('vmax', np.max(field_verts))

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
        colorscale=kwargs['colormap_object'],
        hovertemplate=kwargs['xlabel']+': %{x:.2e}<br>'+\
                        kwargs['ylabel']+': %{y:.2e}<br>'+\
                        kwargs['zlabel']+': %{z:.2e}<br>'+\
                        'field: %{intensity:.2e}',
        name='',
        showscale=False,
        scene=ax['sceneN']
    ))

    if kwargs['colorbar'] and not(ax['colorbar']) and not(kwargs['field_is_nan']):
        ax['colormap_object'] = kwargs['colormap_object']
        fig.add_trace(tool_plotly_colorbar(ax, type='normal'))
        ax['colorbar'] = True

    kwargs['fig'] = fig
    kwargs['ax'] = ax
    tool_set_plot_axis_properties_plotly(self, **kwargs)
    return fig, ax
