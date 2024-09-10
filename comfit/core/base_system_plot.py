from typing import Union, Optional

import numpy as np
import scipy as sp

from comfit.tool import tool_colormap_angle, tool_colormap_bluewhitered, tool_colormap_sunburst
from comfit.tool import tool_create_orthonormal_triad
from comfit.tool import tool_set_plot_axis_properties_matplotlib
from comfit.tool import tool_complete_field

from comfit.plot import plot_surface_matplotlib
from comfit.plot import plot_field_matplotlib
from comfit.plot import plot_field_plotly
from comfit.plot import plot_complex_field_matplotlib
from comfit.plot import plot_complex_field_plotly
from comfit.plot import plot_vector_field_matplotlib
from comfit.plot import plot_vector_field_plotly

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.tri as mtri
import matplotlib.cm as cm
import matplotlib.axes
import matplotlib.figure

import mpl_toolkits

import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.colors as pc

from skimage.measure import marching_cubes

import time

class BaseSystemPlot:
    """ Plotting methods for the base system class"""

    def show(self):
        """Show the plot."""
        if self.plot_lib == 'matplotlib':   
            plt.show()
        elif self.plot_lib == 'plotly':
            go.show() #or fig.show?

    def plot_tool_plotly_add_angle_colorbar3D(self, fig: go.Figure) -> go.Figure:
        """Adds a colorbar to a 3D plot with the angle colormap.

        Args:
            fig: The figure to which the colorbar is added.
        
        Returns:
            The figure with the colorbar added.
        """
        # Generate the custom colormap
        custom_colormap = tool_colormap_angle(pyplot=True)

        fig.add_trace(go.Isosurface(
            x=[self.xmin/self.a0,self.xmax/self.a0], 
            y=[self.ymin/self.a0,self.ymax/self.a0], 
            z=[self.zmin/self.a0, self.zmax/self.a0],
            value=[0,0],
            cmin = -np.pi,
            cmax = np.pi,
            opacity=0,
            colorscale=custom_colormap,
            colorbar=dict(
                tickvals=[-np.pi, -2*np.pi/3, -np.pi/3, 0, np.pi/3, 2*np.pi/3, np.pi],
                ticktext=['-π', '-2π/3', '-π/3', '0', 'π/3', '2π/3', 'π']
            ),
            showscale=True
        ))

        return fig

    def plot_field(self, field, **kwargs):
        return plot_field_plotly(self, field, **kwargs)

    def plot_complex_field(self, complex_field, **kwargs):
        return plot_complex_field_plotly(self, complex_field, **kwargs)

    def plot_angle_field(self, angle_field, **kwargs):
        print('\033[93mWarning: The plot_angle_field function is not yet supported with ploty.\033[0m')

    def plot_vector_field(self, vector_field, **kwargs):
        return plot_vector_field_plotly(self, vector_field, **kwargs)

    def plot_field_in_plane(
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

        # Extend the field if not a complete array is given
        field = tool_complete_field(self, field)

        # Check if the vector field is complex
        if np.iscomplexobj(field):
            print("\033[91mWarning: the provided field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
            print('Max imaginary part: ', np.max(np.imag(field)))
            field = np.real(field)

        # Check if an axis object is provided
        fig = kwargs.get('fig', plt.gcf())
        ax = kwargs.get('ax', None)

        # Kewyord arguments
        colorbar = kwargs.get('colorbar', True)

        # Default values of position and normal vector
        if position is None:
            position = self.rmid

        if normal_vector is None:
            normal_vector=[0,1,0]

        if ax is None:
            ax = fig.add_subplot(111, projection='3d')

        colormap = kwargs.get('colormap', 'viridis')

        if colormap == 'angle':
            colormap = tool_colormap_angle()
        elif colormap == 'bluewhitered':
            colormap = tool_colormap_bluewhitered()
        else:
            colormap = plt.get_cmap(colormap)

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
        colors = colormap(field_normalized)
    
        ax.plot_trisurf((self.xmin+verts[:, 0]*self.dx)/self.a0,
                        (self.ymin+verts[:, 1]*self.dy)/self.a0,
                        faces,
                        (self.zmin+verts[:, 2]*self.dz)/self.a0,
                        facecolor=colors, antialiased=False)

        if colorbar:
            sm = plt.cm.ScalarMappable(cmap=colormap)  
            sm.set_clim(np.min(field_verts),np.max(field_verts))  
            cbar = plt.colorbar(sm, ax=ax, pad=0.2)

        kwargs['grid'] = kwargs.get('grid', True)
        kwargs['ax'] = ax
        tool_set_plot_axis_properties(self,**kwargs)

        return fig, ax

    def plot_complex_field_in_plane(
            self,
            complex_field: np.ndarray,
            normal_vector: Optional[np.ndarray] = None,
            position: Optional[np.ndarray] = None,
            **kwargs
        ) -> matplotlib.axes.Axes:
        """Plots the complex field in a plane perpendicular to the given normal vector using

        Args:
            complex_field (array-like): The complex field to be plotted.
            normal_vector (array-like, optional): The normal vector of the plane. Default is [0,1,0].
            position (array-like, optional): The position of the plane. Default is the middle of the system.
            **kwargs: Keyword arguments for the plot, see https://comfitlib.com/Plotting/.

        Returns:
            The axes containing the plot of the complex field. (matplotlib.axes.Axes)
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
        plot_tool_set_plot_axis_properties(self, **kwargs)

        return fig, ax

    def plot_angle_field_in_plane(
            self,
            angle_field: np.ndarray,
            normal_vector: Optional[np.ndarray] = None,
            position: Optional[np.ndarray] = None,
            **kwargs
        ) -> matplotlib.axes.Axes:
        """Plots the angle field in a plane.

        Args:
            angle_field (numpy.ndarray): The angle field to be plotted.
            normal_vector (array-like, optional): The normal vector of the plane. Default is [0,1,0].
            position (array-like, optional): The position of the plane. Default is the middle of the system.
            **kwargs: Keyword arguments for the plot, see https://comfitlib.com/Plotting/.

        Returns:
            The axes containing the plot. (matplotlib.axes.Axes)
        """

        # Check if the vector field is complex
        if np.iscomplexobj(angle_field):
            print("\033[91mWarning: the angle vector field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
            print('Max imaginary part: ', np.max(np.imag(angle_field)))
            angle_field = np.real(angle_field)

        # Extend the field if not a complete array is given
        angle_field = tool_complete_field(self, angle_field)

        complex_field = np.exp(1j * angle_field)
    
        return self.plot_complex_field_in_plane(complex_field, normal_vector=normal_vector, position=position, **kwargs)

    
    def plot_vector_field_in_plane(
            self,
            vector_field: tuple[np.ndarray, np.ndarray],
            position: Optional[np.ndarray] = None,
            normal_vector: Optional[np.ndarray] = None,
            spacing: int = 2,
            **kwargs
        ) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
        """
        Plots the vector field in a plane.
        
        Args:
            vector_field (tuple): Tuple containing the x and y components of the vector field.
            position (array-like, optional): The position of the plane. Default is the middle of the system.
            normal_vector (array-like, optional): The normal vector of the plane. Default is [0,1,0].
            spacing (int, optional): The spacing for the quiver plot. Default is 5.
            kwargs: Keyword arguments for the plot, see https://comfitlib.com/Plotting/.

        Returns:
            Tuple consisting of
                - The figure containing the plot. matplotlib.figure.Figure.
                - The axes containing the plot. matplotlib.axes.Axes: 
        """

        # Convert the vector field to a numpy array
        vector_field = np.array(vector_field)

        # Extend the field if not a complete array is given
        # TODO: Make this part more efficient (Vidar 11.03.24)
        vector_field_copy = []
        for n in range(vector_field.shape[0]):
            vector_field_copy.append(tool_complete_field(self, vector_field[n]))
        vector_field = np.array(vector_field_copy)

        # Check if the vector field is complex
        if np.iscomplexobj(vector_field):
            print("\033[91mWarning: the provided vector field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
            print('Max imaginary part: ', np.max(np.imag(vector_field)))
            vector_field = np.real(vector_field)
        

        if self.dim != 3:
            raise Exception("The plot in plane function is only defined for 3D fields.")

        # Check if an axis object is provided
        fig = kwargs.get('fig', plt.gcf())
        ax = kwargs.get('ax', None)

        if ax == None:
            ax = fig.add_subplot(111, projection='3d')

        # Default values of position and normal vector
        if position is None:
            position = self.rmid
        
        if normal_vector is None:
            normal_vector = [0,1,0]
        
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

        ax.quiver(x, y, z, U_verts, V_verts, W_verts, color='blue')

        kwargs['ax'] = ax
        plot_tool_set_plot_axis_properties(self, **kwargs)
        return fig, ax