from typing import TYPE_CHECKING, Any, Tuple

if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy as sp
from mpl_toolkits.mplot3d import Axes3D
from skimage.measure import marching_cubes

# Local application imports
from comfit.tool import (
    tool_complete_field,
    tool_colormap,
    tool_set_plot_axis_properties_matplotlib,
    tool_matplotlib_define_2D_plot_ax,
    tool_matplotlib_define_3D_plot_ax
)
from .plot_surface_matplotlib import plot_surface_matplotlib

def plot_complex_field_matplotlib(
        self: 'BaseSystem',
        complex_field: np.ndarray,
        **kwargs: Any
        ) -> Tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """Plot a complex field using matplotlib.

    Creates a visualization of a complex field using different plotting methods
    depending on the dimensionality of the field.

    Parameters
    ----------
    self : BaseSystem
        A BaseSystem (or derived) instance.
    complex_field : np.ndarray
        The complex field to plot
    kwargs : Any
        Keyword arguments for the plot, see https://comfitlib.com/Plotting/

    Returns
    -------
    Tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        The figure and axes objects containing the plot

    Examples
    --------
    >>> system.plot_complex_field_matplotlib(field, colorbar=True)
    """

    kwargs['colormap'] = kwargs.get('colormap', 'angle') # Override the default colormap with 'angle'
    complex_field, fig, ax, kwargs = self.plot_prepare(complex_field, field_type = 'complex', **kwargs)

    # Calculate the magnitude and phase of the complex field
    rho = np.abs(complex_field)
    theta = np.angle(complex_field)

    ###############################################################
    ###################### DIMENSION: 1 ###########################
    ###############################################################

    if self.dim == 1:

        ax = tool_matplotlib_define_2D_plot_ax(fig, ax)

        x = kwargs.get('x', self.x).flatten()

        # Padding for the colorbar.
        padding=0.05

        ax.plot(x, rho, color='black')

        # Color in the graph based on the argument of the complex field
        blend_factor=0.3 # The degree to which the color is blended with white
        cmap = kwargs['colormap_object']

        # Extract coordinates
        dx = x[1] - x[0]
        xmax = x[-1]+dx


        ax.fill_between([x[0],(x[0]+dx/2)], [rho[0],(rho[0]+rho[1])/2],
                        color=(1-blend_factor)*np.array(cmap((theta[0] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1]), 
                        alpha=1)

        for i in range(1,self.xRes-1):
            ax.fill_between([(x[i]-dx/2),x[i]], [(rho[i]+rho[i-1])/2,rho[i]],
                            color=(1-blend_factor)*np.array(cmap((theta[i] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1]), 
                            alpha=1)
            ax.fill_between([x[i],(x[i]+dx/2)], [rho[i],(rho[i]+rho[i+1])/2],
                color=(1-blend_factor)*np.array(cmap((theta[i] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1]),  
                alpha=1)

        ax.fill_between([(xmax-1.5*dx),(xmax-dx)], [(rho[-1]+rho[-2])/2,rho[-1]],
                        color=(1-blend_factor)*np.array(cmap((theta[-1] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1]),  
                        alpha=1)
            
    ###############################################################
    ###################### DIMENSION: 2 ###########################
    ###############################################################

    elif self.dim == 2:
        # Keyword arguments particular to the 2D case
        plot_method = kwargs.get('plot_method', 'phase_angle')

        # Extract coordinates
        x = kwargs.get('x', self.x/self.a0).flatten()
        dx = x[1] - x[0]
        xmax = x[-1]+dx

        y = kwargs.get('y', self.y/self.a0).flatten()
        dy = y[1] - y[0]
        ymax = y[-1]+dy

        # Create a meshgrid
        X, Y = np.meshgrid(x, y, indexing='ij')

        if plot_method == '3Dsurface':
        
            # Padding for the colorbar
            padding=0.2

            # Keyword arguments particular to the 3D surface plot
            grid = kwargs.get('grid', True)
            kwargs['axis_equal'] = False
            
            ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
            kwargs['plot_is_3D'] = True
            
            # Get the colors from a colormap (e.g., hsv, but you can choose any other)
            colors = kwargs['colormap_object']((theta + np.pi) / (2 * np.pi))  # Normalizing theta to [0, 1]

            surf = ax.plot_surface(X, Y, rho, facecolors=colors)

        elif plot_method == 'phase_angle':
            

            # Keyword arguments particular to the phase angle plot
            grid = kwargs.get('grid', False)

            rho_normalized = rho / np.max(rho)

            custom_colormap = kwargs['colormap_object']
            
            # Check if an axis object is provided
            ax = tool_matplotlib_define_2D_plot_ax(fig, ax)

            # Padding for the colorbar
            padding=0.05
            
            mesh = ax.pcolormesh(X, Y, theta, shading='auto', cmap=custom_colormap, vmin=-np.pi, vmax=np.pi)
            mesh.set_alpha(rho_normalized)
 

    elif self.dim == 3:

        grid = kwargs.get('grid', True)
        axis_equal = kwargs.get('axis_equal',True)

        plot_method = kwargs.get('plot_method', 'phase_blob')

        rho_normalized = rho / np.max(rho)

        colormap = kwargs['colormap_object']

        ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
        kwargs['plot_is_3D'] = True
    
        if plot_method == 'phase_angle':
            
            # Padding for the colorbar
            padding=0.2

            for angle in [-2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3]:
                field_to_plot = theta.copy()
                field_to_plot[theta < angle - 1] = float('nan')
                field_to_plot[theta > angle + 1] = float('nan')
                field_to_plot[rho_normalized < 0.01] = float('nan')

                if np.nanmin(field_to_plot) < angle < np.nanmax(field_to_plot):
                    #TODO: make alpha a keyword argument (Vidar 11.03.24)
                    plot_surface_matplotlib(self, field=field_to_plot, 
                                        value=angle, 
                                        color=colormap((angle + np.pi) / (2 * np.pi)), 
                                        alpha=0.5,
                                        ax=ax)

            theta = np.mod(theta, 2 * np.pi)

            field_to_plot = theta.copy()
            field_to_plot[theta < np.pi - 1] = float('nan')
            field_to_plot[theta > np.pi + 1] = float('nan')
            field_to_plot[rho_normalized < 0.01] = float('nan')

            if np.nanmin(field_to_plot) < np.pi < np.nanmax(field_to_plot):
                #TODO: make alpha a keyword argument (Vidar 11.03.24)
                plot_surface_matplotlib(self, field=field_to_plot, 
                                        value=np.pi, 
                                        color=colormap(0), 
                                        alpha=0.5,
                                        ax=ax)
        
        elif plot_method == 'phase_blob':
            # TODO: change the function so that it used plot_surface_matplotlib (Vidar 11.03.24)
            # Padding for the colorbar
            padding=0.2
            
            phase_blob_threshold = kwargs.get('phase_blob_threshold', 0.5)

            if np.nanmin(rho_normalized)<phase_blob_threshold<np.nanmax(rho_normalized):
                verts, faces, _, _ = marching_cubes(rho_normalized, phase_blob_threshold)

                # Calculate the centroids of each triangle
                centroids = np.mean(verts[faces], axis=1)

                # Assuming theta is defined on the same grid as rho
                x, y, z = np.mgrid[0:rho_normalized.shape[0], 0:rho_normalized.shape[1], 0:rho_normalized.shape[2]]

                # Flatten the grid for interpolation
                points = np.c_[x.ravel(), y.ravel(), z.ravel()]
                theta_values = theta.ravel()

                # Interpolate theta at the vertices positions
                theta_faces = sp.interpolate.griddata(points, theta_values, centroids, method='nearest')

                # Normalize theta values for color mapping
                theta_faces_normalized = (theta_faces + np.pi) / (2*np.pi)

                # Map normalized theta values to colors
                colors = colormap(theta_faces_normalized)
                
                # Extract coordinates
                x = kwargs.get('x', self.x/self.a0).flatten()
                dx = x[1] - x[0]
                xmax = x[-1]+dx

                y = kwargs.get('y', self.y/self.a0).flatten()
                dy = y[1] - y[0]
                ymax = y[-1]+dy

                z = kwargs.get('z', self.z/self.a0).flatten()
                dz = z[1] - z[0]
                zmax = z[-1]+dz

                # Plot the complex field
                ax.plot_trisurf(x[0]+verts[:, 0]*dx, 
                                y[0]+verts[:, 1]*dy, 
                                faces, 
                                z[0]+verts[:, 2]*dz, 
                                facecolor=colors, antialiased=False)
            
                # Plot the shadows on the edges
                plot_shadows = kwargs.get('plot_shadows', True)
                if plot_shadows:
                    ax.plot_trisurf((x[0]+0*verts[:, 0]*dx), 
                                    (y[0]+verts[:, 1]*dy), 
                                    faces, 
                                    (z[0]+verts[:, 2]*dz), 
                                    facecolor='black', antialiased=True,
                                    alpha=0.1)

                    ax.plot_trisurf((x[0]+verts[:, 0]*dx), 
                                    (ymax+0*verts[:, 1]*dy), 
                                    faces, 
                                    (z[0]+verts[:, 2]*dz), 
                                    facecolor='black', antialiased=True,
                                    alpha=0.1)
                    
                    ax.plot_trisurf((x[0]+verts[:, 0]*dx), 
                                    (y[0]+verts[:, 1]*dy), 
                                    faces, 
                                    (z[0]+0*verts[:, 2]*dz), 
                                    facecolor='black', antialiased=True,
                                    alpha=0.1)



    if kwargs['colorbar']:

        mappable = plt.cm.ScalarMappable(cmap=kwargs['colormap_object'])
        mappable.set_array([])
        mappable.set_clim(-np.pi, np.pi)
        cbar = plt.colorbar(mappable, ax=ax, pad=padding)

        cticks = kwargs.get('cticks', [-np.pi, -2*np.pi/3, -np.pi/3, 0, np.pi/3, 2*np.pi/3, np.pi])
        cbar.set_ticks(cticks)

        cticklabelse = kwargs.get('cticklabels', [r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$'])
        cbar.set_ticklabels(cticklabelse)

    kwargs['fig'] = fig
    kwargs['ax'] = ax
    tool_set_plot_axis_properties_matplotlib(self, **kwargs)

    return fig, ax
