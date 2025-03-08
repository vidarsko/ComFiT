# Type checking imports
from typing import TYPE_CHECKING, Optional, Any, Tuple

if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

# Standard library imports
import numpy as np
import scipy as sp

# Third-party library imports
from skimage.measure import marching_cubes
import plotly.graph_objects as go

# Local application imports
from comfit.tool import (
    tool_matplotlib_define_3D_plot_ax,
    tool_plotly_define_3D_plot_ax,
    tool_set_plot_axis_properties_matplotlib,
    tool_set_plot_axis_properties_plotly,
    tool_plotly_colorbar
)


def plot_vector_field_in_plane_both_plot_libs(
        self: 'BaseSystem',
        vector_field: np.ndarray,
        position: Optional[np.ndarray] = None,
        normal_vector: Optional[np.ndarray] = None,
        spacing: Optional[int] = None,
        **kwargs: Any
        ) -> Tuple[Any, Any]:
    """Plot a vector field on a plane in 3D space.

    Creates a visualization of a vector field where it intersects with a plane
    defined by a point and a normal vector.

    Parameters
    ----------
    vector_field : np.ndarray
        The vector field to plot. Shape must be (1,xRes,yRes,zRes), 
        (2,xRes,yRes,zRes), or (3,xRes,yRes,zRes).
    position : np.ndarray, optional
        Point through which the plane passes. Defaults to system midpoint.
    normal_vector : np.ndarray, optional
        Normal vector defining the plane orientation. Defaults to [1,1,1].
    spacing : int, optional
        Spacing between vectors in the plot. Defaults to max(xRes//20,1).
    \*\*kwargs : Any
        Additional keyword arguments for plot customization. See https://comfitlib.com/Plotting/.

    Returns
    -------
    Tuple[Any, Any]
        Figure and axes objects from the plotting library.

    Raises
    ------
    Exception
        If the system is not 3D or if the vector field has invalid dimensions.
    """

    plot_lib = kwargs.get('plot_lib', self.plot_lib)

    if self.dim != 3:
        raise Exception("The plot_vector_field_in_plane function is only defined for 3D fields.")

    vector_field, fig, ax, kwargs = self.plot_prepare(vector_field, field_type = 'vector', **kwargs)

    # Default values of position and normal vector
    if position is None:
        position = self.rmid
    
    if normal_vector is None:
        normal_vector = [1,1,1]

    if spacing is None:
        spacing = max(self.xRes//20,1)

    normal_vector = np.array(normal_vector)/np.linalg.norm(normal_vector)
    height_above_plane = (self.x-position[0])*normal_vector[0] + (self.y-position[1])*normal_vector[1] + (self.z-position[2])*normal_vector[2]

    verts, faces, _, _ = marching_cubes(height_above_plane, 0)

    # Assuming field is defined on the same grid as height_above_plane
    x, y, z = np.mgrid[0:height_above_plane.shape[0], 0:height_above_plane.shape[1], 0:height_above_plane.shape[2]]

    # Flatten the grid for interpolation
    points = np.c_[x.ravel(), y.ravel(), z.ravel()]

    # Extract the vector field
    if vector_field.shape == (1,self.xRes,self.yRes,self.zRes):
        n = 1
    elif vector_field.shape == (2,self.xRes,self.yRes,self.zRes):
        n = 2
    elif vector_field.shape == (3,self.xRes,self.yRes,self.zRes):
        n = 3
    else:
        raise Exception("You have entered an invalid field to the plot_vector_field function.")

    U = vector_field[0]
    U_values = U.ravel()
    U_verts = sp.interpolate.griddata(points, U_values, verts, method='nearest')

    if n > 1:
        V = vector_field[1]
        V_values = V.ravel()
        V_verts = sp.interpolate.griddata(points, V_values, verts, method='nearest')
    else:
        V_verts = np.zeros(U_verts.shape)
    
    if n > 2:
        W = vector_field[2]
        W_values = W.ravel()
        W_verts = sp.interpolate.griddata(points, W_values, verts, method='nearest')
    else:
        W_verts = np.zeros(U_verts.shape)

    # Normalize the vectors
    max_vector = np.max(np.sqrt(U_verts ** 2 + V_verts ** 2 + W_verts ** 2))
    U_verts = U_verts / max_vector
    V_verts = V_verts / max_vector
    W_verts = W_verts / max_vector

    # Scale factors
    vx_scale = kwargs.get('vx_scale', spacing*self.dx)
    vy_scale = kwargs.get('vy_scale', spacing*self.dy)
    vz_scale = kwargs.get('vz_scale', spacing*self.dz)

    # Scaling
    U_verts = vx_scale*U_verts
    V_verts = vy_scale*V_verts
    W_verts = vz_scale*W_verts

    x = self.xmin+verts[:, 0]*self.dx
    y = self.ymin+verts[:, 1]*self.dy
    z = self.zmin+verts[:, 2]*self.dz

    # Adjust positions based on the coordinate system
    x = self.xmin + x * self.dx
    y = self.ymin + y * self.dy
    z = self.zmin + z * self.dz

    # Add spacing
    x = x[::spacing]
    y = y[::spacing]
    z = z[::spacing]
    U_verts = U_verts[::spacing]
    V_verts = V_verts[::spacing]
    W_verts = W_verts[::spacing]

    if plot_lib == "plotly":

        ax = tool_plotly_define_3D_plot_ax(fig, ax) #Defines sceneN, plot_dimension

        ax['vmin'] = kwargs.get('vmin', 0)
        ax['vmax'] = kwargs.get('vmax', np.max(max_vector))

        if not kwargs['field_is_nan']:
            fig.add_trace(go.Cone(  x=x.flatten()/self.a0, 
                                    y=x.flatten()/self.a0, 
                                    z=y.flatten()/self.a0, 
                                    u=U_verts.flatten(), 
                                    v=V_verts.flatten(), 
                                    w=W_verts.flatten(), 
                                    colorscale=kwargs['colormap_object'], 
                                    sizemode='scaled', 
                                    sizeref=1, 
                                    showscale=False, 
                                    scene = ax['sceneN'],
                                    hovertemplate='<b>x:</b> %{x:.1f} a₀ <br>' +
                                                '<b>ux:</b> %{u:.1f} <br>' +
                                                '<b>uy:</b> %{v:.1f} <br>' +
                                                '<b>uz:</b> %{w:.1f} <br>' +
                                                '<b>|u|:</b> %{customdata:.1f}<extra></extra>',
                                                customdata=np.sqrt(U_verts**2 + V_verts**2 + W_verts**2).flatten()))

        if kwargs['colorbar'] and not(ax['colorbar']) and not(kwargs['field_is_nan']):
            ax['colormap_object'] = kwargs['colormap_object']
            fig.add_trace(tool_plotly_colorbar(ax, type='normal'))
            ax['colorbar'] = True

        kwargs['fig'] = fig
        kwargs['ax'] = ax
        tool_set_plot_axis_properties_plotly(self, **kwargs)

    elif plot_lib == "matplotlib":
        
        ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
        kwargs['plot_is_3D'] = True

        ax.quiver(x, y, z, U_verts, V_verts, W_verts, color='blue')

        kwargs['fig'] = fig
        kwargs['ax'] = ax
        tool_set_plot_axis_properties_matplotlib(self, **kwargs)
    
    return fig, ax