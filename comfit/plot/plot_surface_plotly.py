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
        alpha: float = 0.5,
        color: str = 'b',
        plt_colormap_object: Any = None,
        hovertemplate: str = None,
        customdata: np.ndarray = None,
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
    alpha : float, optional
        Transparency value. Defaults to 0.5.
    color : str, optional
        Surface color. Defaults to 'b'.
    hovertemplate : str
        Template for the hover text.
    kwargs : Any
        Plotting keyword arguments.

    Returns
    -------
    go.Mesh3D
        The plotly mesh object containing the surface plot.
    """

    customdata = None

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
        
        # Convert colors to 'rgba()' format required by Plotly
        colors = ['rgba({},{},{},{})'.format(*c) for c in (255*colors).astype(int)]
        color = None

    elif isinstance(alpha,np.ndarray):

        alpha_max = np.max(alpha)
        alpha_min = np.min(alpha)

        max_opacity = 0.9

        if (alpha_max - alpha_min)<1e-6:
            alpha = max_opacity*np.ones_like(alpha)
        else:
            alpha = max_opacity*(alpha - alpha_min) / (alpha_max - alpha_min)

        x, y, z = np.mgrid[0:field.shape[0], 0:field.shape[1], 0:field.shape[2]]
        points = np.c_[x.ravel(), y.ravel(), z.ravel()]

        interpolation_method = kwargs.get('interpolation_method', 'linear')
        valid_centroids = ~np.isnan(centroids).any(axis=1)
        rho_normalized_faces = np.full(centroids.shape[0], np.nan)
        rho_normalized_faces[valid_centroids] = sp.interpolate.griddata(points, alpha.ravel(), centroids[valid_centroids], method='nearest')

        colors =   ['rgba({},{},{},{})'.format(*(np.array([
                    (255*color[0]).astype(int),
                    (255*color[1]).astype(int),
                    (255*color[2]).astype(int),
                    0 if np.isnan(rho_normalized_faces[i]) else rho_normalized_faces[i]
                    ]))) for i in range(len(faces))]
        
        hovertemplate=kwargs['xlabel']+': %{x:.2f}<br>'+\
                        kwargs['ylabel']+': %{y:.2f}<br>'+\
                        kwargs['zlabel']+': %{z:.2f}<br>'+\
                        'amplitude: %{customdata:.2e}<br>' + \
                        'phase: '+ f'{value/np.pi:.2f} π '
        
        customdata = alpha_min + rho_normalized_faces*(alpha_max-alpha_min)

        color = 'rgba({},{},{},{})'.format(*(255*np.array(color)).astype(int)) 
        # color = None
        alpha = None

    else:

        colors = None
        color = 'rgba({},{},{},{})'.format(*(255*np.array(color)).astype(int)) 


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


    mesh = go.Mesh3d(
                x=x_new, 
                y=y_new, 
                z=z_new,
                i=faces[:, 0], 
                j=faces[:, 1], 
                k=faces[:, 2],
                facecolor=colors,  # Set color for each face
                color=color,
                showscale=False,
                opacity=alpha,
                hovertemplate=hovertemplate,
                customdata=customdata,
                name='',
                scene=kwargs['ax']['sceneN'],
                hoverlabel=dict(
                    bgcolor=color
                )
            )

    plot_shadows = kwargs.get('plot_shadows', True)
    if plot_shadows:
        shadow_meshes = []

        # Shadow on xy-plane:
        mesh_shadow_xy = go.Mesh3d(
            x=x_new, 
            y=y_new, 
            z=0*z_new + zmin,
            i=faces[:, 0], 
            j=faces[:, 1], 
            k=faces[:, 2],
            color='black',
            showscale=False,
            opacity=0.2,
            name='',
            showlegend=False,
            scene=kwargs['ax']['sceneN']
            )

        # Shadow on xz-plane:
        mesh_shadow_xz = go.Mesh3d(
            x=x_new, 
            y=0*y_new + ymin, 
            z=z_new,
            i=faces[:, 0], 
            j=faces[:, 1], 
            k=faces[:, 2],
            color='black',
            showscale=False,
            opacity=0.2,
            name='',
            showlegend=False,
            scene=kwargs['ax']['sceneN']
            )

        # Shadow on yz-plane:
        mesh_shadow_yz = go.Mesh3d(
            x=0*x_new + xmin, 
            y=y_new, 
            z=z_new,
            i=faces[:, 0], 
            j=faces[:, 1], 
            k=faces[:, 2],
            color='black',
            showscale=False,
            opacity=0.2,
            name='',
            showlegend=False,
            scene=kwargs['ax']['sceneN']
            )

        shadow_meshes.append(mesh_shadow_xy)
        shadow_meshes.append(mesh_shadow_xz)
        shadow_meshes.append(mesh_shadow_yz)
    else:
        shadow_meshes = None

    return mesh, shadow_meshes