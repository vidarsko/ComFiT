import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy as sp
from comfit.tool import tool_complete_field
from comfit.tool import tool_colormap_angle
from comfit.tool import tool_set_plot_axis_properties_matplotlib

from comfit.plot.plot_surface_matplotlib import plot_surface_matplotlib
from mpl_toolkits.mplot3d import Axes3D

from skimage.measure import marching_cubes

def plot_complex_field_matplotlib(self, 
        complex_field: np.ndarray, 
        **kwargs
        ) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """
    Plot a complex field.

    ax=None, plot_method=None, colorbar=False

    Args:
        complex_field (numpy.ndarray): The complex field to plot.
        ax (matplotlib.axes.Axes, optional): The matplotlib axes on which to plot the field.
            If not provided, a new 3D axes will be created.
    
    Returns:
        tuple: A tuple containing:
            - matplotlib.figure.Figure: The figure containing the plot.
            - matplotlib.axes.Axes: The axes containing the plot.
    """

    fig = kwargs.get('fig', plt.gcf())
    ax = kwargs.get('ax', None)

    # Kewyord arguments
    colorbar = kwargs.get('colorbar', True)

    # Extend the field if not a complete array is given
    complex_field = tool_complete_field(self, complex_field)

    # Calculate the magnitude and phase of the complex field
    rho = np.abs(complex_field)
    theta = np.angle(complex_field)

    ###############################################################
    ###################### DIMENSION: 1 ###########################
    ###############################################################

    if self.dim == 1:

        # Keyword arguments particular to the 1D case
        grid = kwargs.get('grid', False)
        axis_equal = kwargs.get('axis_equal',False)
    
        if ax == None:
            fig.clf()
            ax = fig.add_subplot(111)

        # Padding for the colorbar.
        padding=0.05

        ax.plot(self.x/self.a0, rho, color='black')

        # Color in the graph based on the argument of the complex field
        blend_factor=0.3 # The degree to which the color is blended with white
        cmap = tool_colormap_angle()

        ax.fill_between([self.xmin/self.a0,(self.xmin+self.dx/2)/self.a0], [rho[0],(rho[0]+rho[1])/2],
                        color=(1-blend_factor)*np.array(cmap((theta[0] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1]), 
                        alpha=1)

        for i in range(1,self.xRes-1):
            ax.fill_between([(self.x[i]-self.dx/2)/self.a0,self.x[i]/self.a0], [(rho[i]+rho[i-1])/2,rho[i]],
                            color=(1-blend_factor)*np.array(cmap((theta[i] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1]), 
                            alpha=1)
            ax.fill_between([self.x[i]/self.a0,(self.x[i]+self.dx/2)/self.a0], [rho[i],(rho[i]+rho[i+1])/2],
                color=(1-blend_factor)*np.array(cmap((theta[i] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1]),  
                alpha=1)

        ax.fill_between([(self.xmax-1.5*self.dx)/self.a0,(self.xmax-self.dx)/self.a0], [(rho[-1]+rho[-2])/2,rho[-1]],
                        color=(1-blend_factor)*np.array(cmap((theta[-1] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1]),  
                        alpha=1)
            
    ###############################################################
    ###################### DIMENSION: 2 ###########################
    ###############################################################

    elif self.dim == 2:
        # Keyword arguments particular to the 2D case
        plot_method = kwargs.get('plot_method', 'phase_angle')
        axis_equal = kwargs.get('axis_equal',True)

        # Create a meshgrid
        X, Y = np.meshgrid(self.x, self.y, indexing='ij')

        if plot_method == '3Dsurface':
        
            # Padding for the colorbar
            padding=0.2

            # Keyword arguments particular to the 3D surface plot
            grid = kwargs.get('grid', True)
            kwargs['axis_equal'] = False
            
            if ax == None:
                fig.clf()
                ax = fig.add_subplot(111, projection='3d')
            
            # Get the colors from a colormap (e.g., hsv, but you can choose any other)
            colors = tool_colormap_angle()((theta + np.pi) / (2 * np.pi))  # Normalizing theta to [0, 1]

            surf = ax.plot_surface(X/self.a0, Y/self.a0, rho, facecolors=colors)

        elif plot_method == 'phase_angle':
            

            # Keyword arguments particular to the phase angle plot
            grid = kwargs.get('grid', False)

            rho_normalized = rho / np.max(rho)

            custom_colormap = tool_colormap_angle()
            
            # Check if an axis object is provided
            if ax == None:
                fig.clf()
                ax = fig.add_subplot(111)

            # Padding for the colorbar
            padding=0.05
            
            mesh = ax.pcolormesh(X/self.a0, Y/self.a0, theta, shading='auto', cmap=custom_colormap, vmin=-np.pi, vmax=np.pi)
            mesh.set_alpha(rho_normalized)
 

    elif self.dim == 3:

        grid = kwargs.get('grid', True)
        axis_equal = kwargs.get('axis_equal',True)

        plot_method = kwargs.get('plot_method', 'phase_blob')
        
        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        rho_normalized = rho / np.max(rho)

        colormap = tool_colormap_angle()


        ax = kwargs.get('ax', None)

        if ax == None:
            fig.clf()
            ax = fig.add_subplot(111, projection='3d')

        if ax is not None and ax is not isinstance(ax, Axes3D):
            # Get row and column
            subplotspec = ax.get_subplotspec()
            gridspec = subplotspec.get_gridspec()
            row = subplotspec.rowspan.start
            col = subplotspec.colspan.start

            ax.remove()

            fig.add_subplot(gridspec[row, col], projection='3d')
            ax = fig.get_axes()[-1]
    
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
                    
                # Plot the complex field
                ax.plot_trisurf((self.xmin+verts[:, 0]*self.dx)/self.a0, 
                                (self.ymin+verts[:, 1]*self.dy)/self.a0, 
                                faces, 
                                (self.zmin+verts[:, 2]*self.dz)/self.a0, 
                                facecolor=colors, antialiased=False)
            
                # Plot the shadows on the edges
                plot_shadows = kwargs.get('plot_shadows', True)
                if plot_shadows:
                    ax.plot_trisurf((self.xmin+0*verts[:, 0]*self.dx)/self.a0, 
                                    (self.ymin+verts[:, 1]*self.dy)/self.a0, 
                                    faces, 
                                    (self.zmin+verts[:, 2]*self.dz)/self.a0, 
                                    facecolor='black', antialiased=True,
                                    alpha=0.1)

                    ax.plot_trisurf((self.xmin+verts[:, 0]*self.dx)/self.a0, 
                                    (self.ymax+0*verts[:, 1]*self.dy)/self.a0, 
                                    faces, 
                                    (self.zmin+verts[:, 2]*self.dz)/self.a0, 
                                    facecolor='black', antialiased=True,
                                    alpha=0.1)
                    
                    ax.plot_trisurf((self.xmin+verts[:, 0]*self.dx)/self.a0, 
                                    (self.ymin+verts[:, 1]*self.dy)/self.a0, 
                                    faces, 
                                    (self.zmin+0*verts[:, 2]*self.dz)/self.a0, 
                                    facecolor='black', antialiased=True,
                                    alpha=0.1)



    if colorbar:

        mappable = plt.cm.ScalarMappable(cmap=tool_colormap_angle())
        mappable.set_array([])
        mappable.set_clim(-np.pi, np.pi)
        cbar = plt.colorbar(mappable, ax=ax, pad=padding)

        cticks = kwargs.get('cticks', [-np.pi, -2*np.pi/3, -np.pi/3, 0, np.pi/3, 2*np.pi/3, np.pi])
        cbar.set_ticks(cticks)

        cticklabelse = kwargs.get('cticklabels', [r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$'])
        cbar.set_ticklabels(cticklabelse)

    kwargs['ax'] = ax
    kwargs['grid'] = grid
    kwargs['axis_equal'] = axis_equal
    tool_set_plot_axis_properties_matplotlib(self, **kwargs)
    return fig, ax
