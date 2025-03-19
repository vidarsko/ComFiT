from typing import TYPE_CHECKING, Any, Dict

if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem


# Standard library imports
import numpy as np
from skimage.measure import marching_cubes
import scipy as sp

import plotly.graph_objects as go


def plot_surface_plotly(
        self: 'BaseSystem',
        field: np.ndarray,
        value: float,
        ax: Dict,
        alpha: float = 0.5,
        color: str = 'b',
        plt_colormap_object: Any = None,
        **kwargs: Any
        ) -> go.Mesh3d:
    """Plot the surface of the given field.

    Parameters
    ----------
    self : BaseSystem
        A BaseSystem (or derived) instance.
    field : ndarray
        3D array containing the field values
    value : float
        Isosurface value
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. 
    alpha : float, optional
        Transparency value. Defaults to 0.5.
    color : str, optional
        Surface color. Defaults to 'b'.
    curstomdata : np.ndarray
        Custom data to be displayed on hover.
    hovertemplate : str
        Template for the hover text.
    \*\*kwargs : Any
        Plotting keyword arguments.

    Returns
    -------
    go.Mesh3D
        The plotly mesh object containing the surface plot.
    """

    hovertemplate = None

    verts, faces, _, _ = marching_cubes(field, value)
    
    # Calculate the centroids of each triangle
    centroids = np.mean(verts[faces], axis=1)

    # Color according to the angle
    if np.iscomplexobj(color):
        x, y, z = np.mgrid[0:field.shape[0], 0:field.shape[1], 0:field.shape[2]]
        points = np.c_[x.ravel(), y.ravel(), z.ravel()]

        reals = np.real(color)
        imags = np.imag(color)

        interpolation_method = kwargs.get('interpolation_method', 'linear')
        # print("Interpolating points with method ", interpolation_method, ".")
        # print("If this process is slow, consider passing 'interpolation_method='nearest' with the plot_complex_field function.")
        # print("That will speed up the process, but the plot may look less smooth.")
        reals_faces = sp.interpolate.griddata(points, reals.ravel(), centroids, method='nearest')
        imags_faces = sp.interpolate.griddata(points, imags.ravel(), centroids, method='nearest')
        # print("Interpolation done.")

        theta_faces = np.arctan2(imags_faces, reals_faces)

        # Normalize theta values for color mapping
        theta_faces_normalized = (theta_faces + np.pi) / (2*np.pi)

        # Map normalized theta values to colors
        colors = plt_colormap_object(theta_faces_normalized)

        # Convert colors to 'rgba()' format required by Plotly
        color = ['rgba({},{},{},{})'.format(*c) for c in (255*colors).astype(int)]

        theta_faces = np.arctan2(imags_faces, reals_faces)

        hovertemplate=kwargs['xlabel']+': %{x:.2f}<br>'+\
                        kwargs['ylabel']+': %{y:.2f}<br>'+\
                        kwargs['zlabel']+': %{z:.2f}<br>'+\
                        'amplitude: '+ f'{value:.2e}<br>'+\
                        'phase: %{customdata:.2f} π'

        # hovertemplate = kwargs['xlabel']+': %{x:.2f}<br>'+\
        #                 kwargs['ylabel']+': %{y:.2f}<br>'+\
        #                 kwargs['zlabel']+': %{z:.2f}<br>'+\
        #                 'field: '+ f'{value:.2e}'

        customdata = theta_faces/np.pi

    # Extract coordinates
    x = kwargs.get('x', self.x/self.a0).flatten()
    dx = x[1] - x[0]
    xmin = x[0]
    xmax = x[-1]+dx

    y = kwargs.get('y', self.y/self.a0).flatten()
    dy = y[1] - y[0]
    ymin = y[0]
    ymax = y[-1]+dy

    z = kwargs.get('z', self.z/self.a0).flatten()
    dz = z[1] - z[0]
    zmin = z[0]
    zmax = z[-1]+dz

    # Create the mesh object
    x_new = (verts[:, 0] * dx + xmin) 
    y_new = (verts[:, 1] * dy + ymin) 
    z_new = (verts[:, 2] * dz + zmin) 


    if hovertemplate is None:
        hovertemplate = kwargs['xlabel']+': %{x:.2f}<br>'+\
                        kwargs['ylabel']+': %{y:.2f}<br>'+\
                        kwargs['zlabel']+': %{z:.2f}<br>'+\
                        'field: '+ f'{value:.2e}'

    print(hovertemplate)
    mesh = go.Mesh3d(
                x=x_new, 
                y=y_new, 
                z=z_new,
                i=faces[:, 0], 
                j=faces[:, 1], 
                k=faces[:, 2],
                facecolor=color,  # Set color for each face
                showscale=False,
                hovertemplate=hovertemplate,
                customdata=customdata,
                name='',
                scene=ax['sceneN']
            )

    return mesh