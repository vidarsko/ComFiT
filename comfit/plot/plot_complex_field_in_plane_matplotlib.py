import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from typing import Optional
from comfit.tool.tool_complete_field import tool_complete_field
from comfit.tool.tool_set_plot_axis_properties_matplotlib import tool_set_plot_axis_properties_matplotlib

from comfit.tool.tool_colormaps import tool_colormap_angle

from skimage.measure import marching_cubes
import matplotlib

def plot_complex_field_in_plane_matplotlib(
            self,
            complex_field: np.ndarray,
            normal_vector: Optional[np.ndarray] = None,
            position: Optional[np.ndarray] = None,
            **kwargs
        ) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
        """Plots the complex field in a plane perpendicular to the given normal vector using

        Args:
            complex_field (array-like): The complex field to be plotted.
            normal_vector (array-like, optional): The normal vector of the plane. Default is [0,1,0].
            position (array-like, optional): The position of the plane. Default is the middle of the system.
            **kwargs: Keyword arguments for the plot, see https://comfitlib.com/Plotting/.

        Returns:
            tuple: A tuple containing:
                - matplotlib.figure.Figure: The figure containing the plot.
                - matplotlib.axes.Axes: The axes containing the plot
        """

        if self.dim != 3:
            raise Exception("The plot in plane function is only defined for 3D fields.")

        # Extend the field if not a complete array is given
        complex_field = tool_complete_field(self, complex_field)

        # Default values of position and normal vector
        if position is None:
            position = self.rmid

        if normal_vector is None:
            normal_vector = [0,1,0]

        # Calculate the magnitude and phase of the complex field
        rho = np.abs(complex_field)
        theta = np.angle(complex_field)
        
        # Check if an axis object is provided
        fig = kwargs.get('fig', plt.gcf())
        ax = kwargs.get('ax', None)

        if ax == None:
            fig.clf()
            ax = fig.add_subplot(111, projection='3d')

        # Kewyord arguments
        colorbar = kwargs.get('colorbar', True)

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
        colormap = tool_colormap_angle()
        colors = colormap(theta_normalized)

        # Blend the colors with white according to rho (normalized)
        colors[:,3] = (rho_verts/np.max(rho_verts)).ravel()
    
        ax.plot_trisurf((self.xmin+verts[:, 0]*self.dx)/self.a0,
                        (self.ymin+verts[:, 1]*self.dy)/self.a0,
                        faces,
                        (self.zmin+verts[:, 2]*self.dz)/self.a0,
                        facecolor=colors, antialiased=True)

        # Create a colorbar
        if colorbar:
            mappable = plt.cm.ScalarMappable(cmap=tool_colormap_angle())
            mappable.set_array([])
            mappable.set_clim(-np.pi, np.pi)
            cbar = plt.colorbar(mappable, ax=ax, pad=0.2)
            cbar.set_ticks(np.array([-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi]))
            cbar.set_ticklabels([r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$'])

        kwargs['grid'] = kwargs.get('grid', True)
        kwargs['ax'] = ax
        tool_set_plot_axis_properties_matplotlib(self, **kwargs)

        return fig, ax