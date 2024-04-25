import numpy as np
import scipy as sp

from comfit.tools.tool_colormaps import tool_colormap_angle, tool_colormap_bluewhitered, tool_colormap_sunburst
from comfit.tools.tool_create_orthonormal_triad import tool_create_orthonormal_triad

import matplotlib.pyplot as plt
from skimage.measure import marching_cubes
import matplotlib.colors as mcolors
import matplotlib.tri as mtri

class BaseSystemPlot:
    # PLOTTING FUNCTIONS

    def plot_tool_set_axis_properties(self, **kwargs):
        """
        Sets the properties of the axis for a plot.
        
        Args:
            kwargs: keyword arguments for the axis properties
        """

        # Check if an axis object is provided
        if 'ax' in kwargs:
            ax = kwargs['ax']

        # Set the xlabel
        if 'xlabel' in kwargs:
            ax.set_xlabel(kwargs['xlabel'])
        else:
            ax.set_xlabel('$x/a_0$')

        if 'title' in kwargs:
                ax.set_title(kwargs['title'])

        #TODO: Check if I should instead include a fig object somewhere (Vidar 05.02.24)
        if 'suptitle' in kwargs:
            plt.suptitle(kwargs['suptitle'])

        # Set the grid
        if 'grid' in kwargs:
            ax.grid(kwargs['grid'])
        else:
            ax.grid(True)

        xlim = [self.xmin/self.a0, (self.xmax-self.dx)/self.a0]

        if 'xmin' in kwargs:
            xlim[0] = kwargs['xmin'] / self.a0

        if 'xmax' in kwargs:
            xlim[1] = kwargs['xmax'] / self.a0

        if 'xlim' in kwargs:
            xlim = np.array(kwargs['xlim']) / self.a0

        ylim = [self.ymin/self.a0, (self.ymax-self.dy)/self.a0] if self.dim > 1 else None

        if 'ymin' in kwargs:
            ylim[0] = kwargs['ymin'] / self.a0 if self.dim > 1 else kwargs['ymin']
        
        if 'ymax' in kwargs:
            ylim[1] = kwargs['ymax'] / self.a0 if self.dim > 1 else kwargs['ymax']
        
        if 'ylim' in kwargs:
            ylim = np.array(kwargs['ylim'])/self.a0 if self.dim > 1 else kwargs['ylim']

        zlim = [self.zmin/self.a0, (self.zmax-self.dz)/self.a0] if self.dim > 2 else None

        if 'zmin' in kwargs:
             zlim[0] = kwargs['zmin'] / self.a0 if self.dim > 2 else kwargs['zmin']

        if 'zmax' in kwargs:
            zlim[1] = kwargs['zmax'] / self.a0 if self.dim > 2 else kwargs['zmax']

        if 'zlim' in kwargs:
            zlim = np.array(kwargs['zlim'])/self.a0 if self.dim > 2 else kwargs['zlim']

        # Custom x-ticks
        if 'xticks' in kwargs:
            ax.set_xticks(kwargs['xticks'])
        if 'xticklabels' in kwargs:
            ax.set_xticklabels(kwargs['xticklabels'])

        # Custom y-ticks
        if 'yticks' in kwargs:
            ax.set_yticks(kwargs['yticks'])
        if 'yticklabels' in kwargs:
            ax.set_yticklabels(kwargs['yticklabels'])
        
        # Custom z-ticks
        if 'zticks' in kwargs:
            ax.set_zticks(kwargs['zticks'])
        if 'zticklabels' in kwargs:
            ax.set_zticklabels(kwargs['zticklabels'])
        
        if self.dim == 1:

            # Set the ylabel
            if 'ylabel' in kwargs:
                ax.set_ylabel(kwargs['ylabel'])

            # Set the title
            if 'title' in kwargs:
                ax.set_title(kwargs['title'])
            
            ax.set_xlim(xlim[0], xlim[1])

            if ylim is not None:
                ax.set_ylim(ylim[0], ylim[1])

            if zlim is not None:
                ax.set_zlim(zlim[0], zlim[1])
            
        
        elif self.dim == 2:

            if 'ylabel' in kwargs:
                ax.set_ylabel(kwargs['ylabel'])
            else:
                ax.set_ylabel('$y/a_0$')
            
            ax.set_xlim(xlim[0], xlim[1])
            ax.set_ylim(ylim[0], ylim[1])
            if zlim is not None:
                ax.set_zlim(zlim[0], zlim[1])

            # Set the aspect ratio
            axis_equal = kwargs.get('axis_equal', True)
            if axis_equal:
                ax.set_aspect('equal')

        elif self.dim == 3:

            if 'ylabel' in kwargs:
                ax.set_ylabel(kwargs['ylabel'])
            else:
                ax.set_ylabel('$y/a_0$')

            if 'zlabel' in kwargs:
                ax.set_zlabel(kwargs['zlabel'])
            else:
                ax.set_zlabel('$z/a_0$')

            ax.set_xlim3d(xlim[0], xlim[1])
            ax.set_ylim3d(ylim[0], ylim[1])
            ax.set_zlim3d(zlim[0], zlim[1])

            # Set the aspect ratio
            axis_equal = kwargs.get('axis_equal', True)
            if axis_equal:
                ax.set_aspect('equal')

        return ax

    def plot_tool_extend_field(self,field):
        """
        Extends the field in case not a complete array is given. 
        For instance, if a field in 3 dimensions is calculated using only x, then the field is extended to the full 3D array.
        Args:
            field: The field to be extended
        Returns:
            np.ndarray: The extended field
        """    
        # 2 dimensional fields
        if field.shape == (self.xRes,1):
            field = np.tile(field,(1,self.yRes))

        elif field.shape == (1,self.yRes):
            field = np.tile(field,(self.xRes,1))

        # 3 dimensional fields
        elif field.shape == (self.xRes,1,1):
            field = np.tile(field,(1,self.yRes,1))
            field = np.tile(field,(1,1,self.zRes))

        elif field.shape == (1,self.yRes,1):
            field = np.tile(field,(self.xRes,1,1))
            field = np.tile(field,(1,1,self.zRes))

        elif field.shape == (1,1,self.zRes):
            field = np.tile(field,(self.xRes,1,1))
            field = np.tile(field,(1,self.yRes,1))

        elif field.shape == (self.xRes,self.yRes,1):
            field = np.tile(field,(1,1,self.zRes))

        elif field.shape == (self.xRes,1,self.zRes):
            field = np.tile(field,(1,self.yRes,1))

        elif field.shape == (1,self.yRes,self.zRes):
            field = np.tile(field,(self.xRes,1,1))

        return field

    def plot_tool_surface(self, **kwargs):
        """
        Plots the surface of the given field.
        """
        field = kwargs['field']
        value = kwargs['value']
        ax = kwargs.get('ax', plt.gca())
        alpha = kwargs.get('alpha', 0.5)
        color = kwargs.get('color', 'b')

        verts, faces, _, _ = marching_cubes(field, value)
        x = (self.xmin+verts[:, 0]*self.dx)/self.a0
        y = (self.ymin+verts[:, 1]*self.dy)/self.a0
        z = (self.zmin+verts[:, 2]*self.dz)/self.a0
        ax.plot_trisurf(x, y, faces, z, alpha=alpha, color=color)
    
    def plot_field(self, field, **kwargs):
        """
        Plots the given (real) field.
        
        Args:
            field (array-like): The field to be plotted.
            **kwargs: Keyword arguments for the plot.
                See github.com/vidarsko/ComFiT/blob/main/docs/ClassBaseSystem.md 
                for a full list of keyword arguments.
        
        Returns:
            matplotlib.axes.Axes: The axes containing the plot.
        """

        if field.dtype == bool:
            field = field.astype(float)

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

        # Extend the field if not a complete array is given
        field = self.plot_tool_extend_field(field)
        
        if self.dim == 1:

            # Keyword arguments particular to the 1D case
            kwargs['grid'] = kwargs.get('grid', True)

            if ax == None:
                fig.clf()
                ax = fig.add_subplot(111)

            ax.plot(self.x/self.a0, field)


        if self.dim == 2:
            
            # Keyword arguments particular to the 2D case
            kwargs['grid'] = kwargs.get('grid', False)

            if ax == None:
                fig.clf()
                ax = fig.add_subplot(111)
            
            # Set the colormap
            colormap = kwargs.get('colormap', 'viridis')

            if colormap == 'bluewhitered':
                colormap = tool_colormap_bluewhitered()

            elif colormap == 'sunburst':
                colormap = tool_colormap_sunburst()
            else:
                colormap = plt.get_cmap(colormap)

            # Value limits symmetric
            vlim_symmetric = kwargs.get('vlim_symmetric', False)

            X, Y = np.meshgrid(self.x, self.y, indexing='ij')

            pcm = ax.pcolormesh(X / self.a0, Y / self.a0, field, shading='gouraud', cmap=colormap)

            xlim = [self.xmin, self.xmax-self.dx]
            ylim = [self.ymin, self.ymax-self.dy]

            limits_provided = False
            if 'xlim' in kwargs:
                xlim = kwargs['xlim']
                limits_provided = True
            else:
                if 'xmin' in kwargs:
                    xlim[0] = kwargs['xmin']
                    limits_provided = True
                
                if 'xmax' in kwargs:
                    xlim[1] = kwargs['xmax']
                    limits_provided = True

            if 'ylim' in kwargs:
                ylim = kwargs['ylim']
                limits_provided = True
            else:
                if 'ymin' in kwargs:
                    ylim[0] = kwargs['ymin']
                    limits_provided = True
                    
                if 'ymax' in kwargs:
                    ylim[1] = kwargs['ymax']
                    limits_provided = True

            # If explicit limits are provided, use them to change the vlim ranges
            if limits_provided:
                region_to_plot = np.zeros(self.dims).astype(bool)
                region_to_plot[(xlim[0] <= X)*(X <= xlim[1])*(ylim[0] <= Y)*(Y <= ylim[1])] = True
                vlim = [np.min(field[region_to_plot]), np.max(field[region_to_plot])]

            else:
                vlim = [np.min(field), np.max(field)]
            
            # Set the value limitses
            if 'vlim' in kwargs:
                vlim = kwargs['vlim']
            else:
                if 'vmin' in kwargs:
                    vlim[0] = kwargs['vmin']
                if 'vmax' in kwargs:
                    vlim[1] = kwargs['vmax']

            if vlim[1] - vlim[0] < 1e-10:
                vlim = [vlim[0]-0.05, vlim[1]+0.05]

            pcm.set_clim(vmin=vlim[0], vmax=vlim[1])

            if 'vlim_symmetric' in kwargs:
                vlim_symmetric = kwargs['vlim_symmetric']
                if vlim_symmetric:
                    cmax = abs(field).max()
                    cmin = -cmax
                    pcm.set_clim(vmin=cmin, vmax=cmax)

            colorbar = kwargs.get('colorbar', True)

            if colorbar:
                cbar = plt.colorbar(pcm, ax=ax)

        elif self.dim == 3:

            if ax == None:
                plt.clf()
                ax = plt.gcf().add_subplot(111, projection='3d')

            field_min = np.min(field)
            field_max = np.max(field)

            number_of_layers = kwargs.get('number_of_layers', 1)
            alpha = kwargs.get('alpha', 0.5)
            
            if 'vlim' in kwargs:
                vlim = kwargs['vlim']
                vmin = vlim[0]
                vmax = vlim[1]
            else:
                vmin = field_min
                vmax = field_max

            if 'layer_values' in kwargs:
                layer_values = np.concatenate([[-np.inf], kwargs['layer_values'], [np.inf]])
            else: 
                layer_values = np.linspace(vmin, vmax, number_of_layers + 2)

            if 'colormap' in kwargs:
                colormap = kwargs['colormap']
                if colormap == 'bluewhitered':
                    colormap = tool_colormap_bluewhitered()

                elif colormap == 'sunburst':
                    colormap = tool_colormap_sunburst()

                else:
                    colormap = plt.get_cmap(colormap)
            else: 
                colormap = plt.get_cmap('viridis')
            

            if field_min < layer_values[1] < field_max:
                self.plot_tool_surface(field=field, 
                                        value=layer_values[1], 
                                        color=colormap((layer_values[1]-vmin) / (vmax-vmin)), 
                                        alpha=alpha,
                                        ax=ax)

            for layer_value in layer_values[2:-1]:
                if field_min < layer_value < field_max:
                    self.plot_tool_surface(field=field, 
                                            value=layer_value, 
                                            color=colormap((layer_value-vmin) / (vmax-vmin)), 
                                            alpha=alpha,
                                            ax=ax)


            if 'colorbar' in kwargs:
                colorbar = kwargs['colorbar']
            else:
                colorbar = True

            if colorbar:
                sm = plt.cm.ScalarMappable(cmap=colormap)
                sm.set_clim(vmin, vmax)
                plt.colorbar(sm, ax=ax, pad=0.2)
        
        kwargs['ax'] = ax
        self.plot_tool_set_axis_properties(**kwargs)
        return fig, ax

    def plot_complex_field(self, complex_field, **kwargs):
        """
        Plot a complex field.

        ax=None, plot_method=None, colorbar=False

        Args:
            complex_field (numpy.ndarray): The complex field to plot.
            ax (matplotlib.axes.Axes, optional): The matplotlib axes on which to plot the field.
                If not provided, a new 3D axes will be created.
        
        Returns:
            matplotlib.axes.Axes: The axes containing the plot.

        """

        # Check if an axis object is provided
        fig = kwargs.get('fig', plt.gcf())
        ax = kwargs.get('ax', None)

        # Kewyord arguments
        colorbar = kwargs.get('colorbar', True)

        # Extend the field if not a complete array is given
        complex_field = self.plot_tool_extend_field(complex_field)

        # Calculate the magnitude and phase of the complex field
        rho = np.abs(complex_field)
        theta = np.angle(complex_field)

        if self.dim == 1:

            # Keyword arguments particular to the 1D case
            grid = kwargs.get('grid', False)

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


        elif self.dim == 2:
            # Keyword arguments particular to the 2D case
            plot_method = kwargs.get('plot_method', 'phase_angle')

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
                
                # Padding for the colorbar
                padding=0.05

                # Keyword arguments particular to the phase angle plot
                grid = kwargs.get('grid', False)
                
                # Check if an axis object is provided
                if ax == None:
                    fig.clf()
                    ax = fig.add_subplot(111)

                rho_normalized = rho / np.max(rho)
                custom_colormap = tool_colormap_angle()

                mesh = ax.pcolormesh(X/self.a0, Y/self.a0, theta, shading='auto', cmap=custom_colormap, vmin=-np.pi, vmax=np.pi)
                mesh.set_alpha(rho_normalized)

        elif self.dim == 3:

            grid = kwargs.get('grid', True)

            plot_method = kwargs.get('plot_method', 'phase_blob')
            
            X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

            rho_normalized = rho / np.max(rho)

            colormap = tool_colormap_angle()

            ax = kwargs.get('ax', None)

            if ax == None:
                fig.clf()
                ax = fig.add_subplot(111, projection='3d')
        
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
                        self.plot_tool_surface(field=field_to_plot, 
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
                    self.plot_tool_surface(field=field_to_plot, 
                                            value=np.pi, 
                                            color=colormap(0), 
                                            alpha=0.5,
                                            ax=ax)
            
            elif plot_method == 'phase_blob':
                # TODO: change the function so that it used plot_tool_surface (Vidar 11.03.24)
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


        # Create a colorbar
        if colorbar:
            mappable = plt.cm.ScalarMappable(cmap=tool_colormap_angle())
            mappable.set_array([])
            mappable.set_clim(-np.pi, np.pi)
            cbar = plt.colorbar(mappable, ax=ax, pad=padding)
            cbar.set_ticks(np.array([-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi]))
            cbar.set_ticklabels([r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$'])

        kwargs['ax'] = ax
        kwargs['grid'] = grid
        self.plot_tool_set_axis_properties(**kwargs)
        return fig, ax

    def plot_angle_field(self, angle_field, **kwargs):
        """
        Plot the angle field.

        Args:
            field (array-like): The angle field values.
            ax (matplotlib.axes.Axes, optional): The axes to plot the angle field on. If not provided, a new subplot will be created.
        
        Returns:
            matplotlib.axes.Axes: The axes containing the plot.
        """

        # Check if the vector field is complex
        if np.iscomplexobj(angle_field):
            print("\033[91mWarning: the angle vector field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
            print('Max imaginary part: ', np.max(np.imag(angle_field)))
            angle_field = np.real(angle_field)

        # Extend the field if not a complete array is given
        angle_field = self.plot_tool_extend_field(angle_field)

        # Normalize around 0
        angle_field = np.mod(angle_field, 2 * np.pi) - np.pi        

        if self.dim == 1:
            if 'vlim' in kwargs:
                vlim = kwargs['vlim']
            else:
                kwargs['vlim'] = [-np.pi, np.pi]
                kwargs['yticks'] = [-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi]
                kwargs['yticklabels'] = [r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$']

            
            
            return self.plot_field(angle_field, **kwargs)
        
        elif self.dim > 1:
            complex_field = np.exp(1j * angle_field)

            kwargs['plot_method'] = 'phase_angle'

            return self.plot_complex_field(complex_field, **kwargs)

    def plot_vector_field(self, vector_field, spacing=5, **kwargs):
        """
        Plots a vector field.

        Args:
        vector_field (tuple): Tuple containing the x and y components of the vector field.
        spacing (int, optional): The spacing for the quiver plot. Default is 5.
        kwargs: Keyword arguments for the plot, see https://vidarsko.github.io/ComFiT/Plotting/.

        Returns:
        matplotlib.figure.Figure: The figure containing the plot.
        matplotlib.axes.Axes: The axes containing the plot.
        """
        # Convert the vector field to a numpy array
        vector_field = np.array(vector_field)
        
        # Extend the field if not a complete array is given
        # TODO: Make this part more efficient (Vidar 11.03.24)
        vector_field_copy = []
        for n in range(vector_field.shape[0]):
            vector_field_copy.append(self.plot_tool_extend_field(vector_field[n]))
        vector_field = np.array(vector_field_copy)

        # Check if the vector field is complex
        if np.iscomplexobj(vector_field):
            print("\033[91mWarning: the provided vector field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
            print('Max imaginary part: ', np.max(np.imag(vector_field)))
            vector_field = np.real(vector_field)

        # Check if an axis object is provided
        fig = kwargs.get('fig', plt.gcf())
        ax = kwargs.get('ax', None)

        def add_spacing_2D(X,Y,U,V,spacing):
            X = X[::spacing, ::spacing]
            Y = Y[::spacing, ::spacing]
            U = U[::spacing, ::spacing]
            V = V[::spacing, ::spacing]
            return X,Y,U,V

        def add_spacing_3D(X,Y,Z,U,V,W,spacing):
            X = X[::spacing, ::spacing, ::spacing]
            Y = Y[::spacing, ::spacing, ::spacing]
            Z = Z[::spacing, ::spacing, ::spacing]
            U = U[::spacing, ::spacing, ::spacing]
            V = V[::spacing, ::spacing, ::spacing]
            W = W[::spacing, ::spacing, ::spacing]
            return X,Y,Z,U,V,W            

        if self.dim == 1:
            
            if vector_field.shape == (1,self.xRes):
                if ax == None:
                    ax = plt.gcf().add_subplot(111)

                X, Y = np.meshgrid(self.x, np.array([0]), indexing='ij')

                U = np.zeros(X.shape)
                V = np.zeros(X.shape)

                V[:,0] = vector_field[0]

                X,Y,U,V = add_spacing_2D(X,Y,U,V,spacing)

                ax.quiver(X/self.a0, Y/self.a0, U, V, color='blue', angles='xy', scale_units='xy', scale=1)
                
                kwargs['ylim'] = [np.min(vector_field[0]), np.max(vector_field[0])]

            elif vector_field.shape == (2,self.xRes):
                if ax == None:
                    fig.clf()
                    ax = fig.add_subplot(111, projection='3d')
                
                X, Y, Z = np.meshgrid(self.x, np.array([0]), np.array([0]), indexing='ij')

                U = np.zeros(X.shape)
                V = np.zeros(X.shape)
                W = np.zeros(X.shape)

                V[:,0,0] = vector_field[0]
                W[:,0,0] = vector_field[1]

                X,Y,Z,U,V,W = add_spacing_3D(X,Y,Z,U,V,W,spacing)

                ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')

                kwargs['ylim'] = [np.min(vector_field[0]), np.max(vector_field[0])]
                delta_y = kwargs['ylim'][1] - kwargs['ylim'][0]

                kwargs['zlim'] = [np.min(vector_field[1]), np.max(vector_field[1])]
                delta_z = kwargs['zlim'][1] - kwargs['zlim'][0]

                if delta_y < 0.15*delta_z:
                    kwargs['ylim'] = [kwargs['ylim'][0] - 0.15*delta_z, kwargs['ylim'][1] + 0.15*delta_z]
                
                if delta_z < 0.15*delta_y:
                    kwargs['zlim'] = [kwargs['zlim'][0] - 0.15*delta_y, kwargs['zlim'][1] + 0.15*delta_y]

            elif vector_field.shape == (3,self.xRes):

                if ax == None:
                    fig.clf()
                    ax = fig.add_subplot(111, projection='3d')
                
                X, Y, Z = np.meshgrid(self.x, np.array([0]), np.array([0]), indexing='ij')

                U = np.zeros(X.shape)
                V = np.zeros(X.shape)
                W = np.zeros(X.shape)

                # Add the vector field
                U[:,0,0] = vector_field[0]
                V[:,0,0] = vector_field[1]
                W[:,0,0] = vector_field[2]

                # Add spacing
                X,Y,Z,U,V,W = add_spacing_3D(X,Y,Z,U,V,W,spacing)

                # Normalize the vectors
                max_vector = np.max(np.sqrt(U ** 2 + V ** 2 + W ** 2))

                # Normalizing
                U = U/max_vector
                V = V/max_vector
                W = W/max_vector

                # Scale factors
                vx_scale = kwargs.get('vx_scale', spacing*self.xmax/15/self.a0)
                vy_scale = kwargs.get('vy_scale', 1)
                vz_scale = kwargs.get('vz_scale', 1)

                # Scaling
                U = vx_scale*U
                V = vy_scale*V
                W = vz_scale*W

                ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')
                
                kwargs['ylim'] = [-1,1]
                kwargs['zlim'] = [-1,1]

            else:
                raise Exception("You have entered an invalid field to the plot_vector_field function.")

        elif self.dim == 2:

            X, Y = np.meshgrid(self.x, self.y, indexing='ij')

            if vector_field.shape == (1,self.xRes,self.yRes):

                if ax == None:
                    ax = plt.gcf().add_subplot(111)

                X, Y = np.meshgrid(self.x, self.y, indexing='ij')

                U = vector_field[0]
                V = np.zeros(X.shape)

                X,Y,U,V = add_spacing_2D(X,Y,U,V,spacing)

                # Normalize the vectors
                max_vector = np.max(np.sqrt(U ** 2))
                
                U = U / max_vector

                # Scale factors
                vx_scale = kwargs.get('vx_scale', spacing*self.dx/self.a0)

                # Scaling
                U = vx_scale*U

                ax.quiver(X/self.a0, Y/self.a0, U, V, color='blue',
                    angles='xy', scale_units='xy', scale=1,
                    headwidth=3, headlength=4, headaxislength=3, pivot='middle')

            elif vector_field.shape == (2,self.xRes,self.yRes):

                if ax == None:
                    ax = plt.gcf().add_subplot(111)

                X, Y, U, V = add_spacing_2D(X,Y,vector_field[0],vector_field[1],spacing)


                max_vector = np.max(np.sqrt(U ** 2 + V ** 2))
                
                U = U / max_vector
                V = V / max_vector

                # Scale factors
                vx_scale = kwargs.get('vx_scale', 1.0*spacing*self.dx/self.a0)
                vy_scale = kwargs.get('vy_scale', 1.0*spacing*self.dy/self.a0)

                # Scaling
                U = vx_scale*U
                V = vy_scale*V

                ax.quiver(X/self.a0, Y/self.a0, U, V, color='blue', 
                    angles='xy', scale_units='xy', scale=1,
                    headwidth=3, headlength=4, headaxislength=3, pivot='middle')

            elif vector_field.shape == (3,self.xRes,self.yRes):

                if ax == None:
                    fig.clf()
                    ax = fig.add_subplot(111, projection='3d')

                X, Y, Z = np.meshgrid(self.x, self.y, np.array([0]), indexing='ij')
                U = np.zeros(X.shape)
                V = np.zeros(X.shape)
                W = np.zeros(X.shape)
                
                U[:,:,0] = vector_field[0]
                V[:,:,0] = vector_field[1]
                W[:,:,0] = vector_field[2]

                X,Y,Z,U,V,W = add_spacing_3D(X,Y,Z,U,V,W,spacing)

                max_vector = np.max(np.sqrt(U ** 2 + V ** 2 + W ** 2))

                # Normalizing
                U = U / max_vector
                V = V / max_vector
                W = W / max_vector
                
                # Scale factors
                vx_scale = kwargs.get('vx_scale', 2*spacing*self.xmax/max_vector/self.a0)
                vy_scale = kwargs.get('vy_scale', 2*spacing*self.ymax/max_vector/self.a0)
                vz_scale = kwargs.get('vz_scale', spacing/self.a0)

                # Scaling
                U = vx_scale*U
                V = vy_scale*V
                W = vz_scale*W

                ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')

                kwargs['axis_equal'] = False
                kwargs['zlim'] = [-spacing,spacing]

        elif self.dim == 3:

            X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

            if ax == None:
                fig.clf()
                ax = fig.add_subplot(111, projection='3d')

            if vector_field.shape == (1,self.xRes,self.yRes,self.zRes):
                # Define the vector field
                U = vector_field[0]
                V = np.zeros(U.shape)
                W = np.zeros(U.shape)

                # Add spacing
                X,Y,Z,U,V,W = add_spacing_3D(X,Y,Z,U,V,W,spacing)

                # Normalize the vectors
                max_vector = np.max(np.sqrt(U ** 2))
                U = U / max_vector

                # Scale factors
                vx_scale = kwargs.get('vx_scale', spacing*self.dx/self.a0)

                # Scaling
                U = vx_scale*U

                ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')

            elif vector_field.shape == (2,self.xRes,self.yRes,self.zRes):
                
                # Define the vector field
                U = vector_field[0]
                V = vector_field[1]
                W = np.zeros(U.shape)

                # Add spacing
                X,Y,Z,U,V,W = add_spacing_3D(X,Y,Z,U,V,W,spacing)

                # Normalize the vectors
                max_vector = np.max(np.sqrt(U ** 2 + V ** 2))
                U = U / max_vector
                V = V / max_vector

                # Scale factors
                vx_scale = kwargs.get('vx_scale', spacing*self.dx/self.a0)
                vy_scale = kwargs.get('vy_scale', spacing*self.dy/self.a0)

                # Scaling
                U = vx_scale*U
                V = vy_scale*V

                ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')

            elif vector_field.shape == (3,self.xRes,self.yRes,self.zRes):
                
                # Define the vector field
                U = vector_field[0]
                V = vector_field[1]
                W = vector_field[2]

                X,Y,Z,U,V,W = add_spacing_3D(X,Y,Z,U,V,W,spacing)

                # Normalize the vectors
                max_vector = np.max(np.sqrt(U ** 2 + V ** 2 + W ** 2))
                U = U / max_vector
                V = V / max_vector
                W = W / max_vector

                # Scale factors
                vx_scale = kwargs.get('vx_scale', spacing*self.dx/self.a0)
                vy_scale = kwargs.get('vy_scale', spacing*self.dy/self.a0)
                vz_scale = kwargs.get('vz_scale', spacing*self.dz/self.a0)

                # Scaling
                U = vx_scale*U
                V = vy_scale*V
                W = vz_scale*W

                ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')

        kwargs['ax'] = ax
        self.plot_tool_set_axis_properties(**kwargs)
        return fig, ax

    def plot_field_in_plane(self, field, normal_vector=None, position=None, 
                        **kwargs):
        """
        Plots the field in a plane perpendicular to the given normal vector using
        scipy.interpolate.griddata and plt.plot_trisurf.

        Args:
            field (array-like): The field to be plotted.
            normal_vector (array-like, optional): The normal vector of the plane. Default is [0,1,0].
            position (array-like, optional): The position of the plane. Default is the middle of the system.
            **kwargs: Keyword arguments for the plot, see https://vidarsko.github.io/ComFiT/Plotting/. 
        
        Returns:
            matplotlib.axes.Axes: The axes containing the plot.
        """

        if self.dim != 3:
            raise Exception("The plot in plane function is only defined for 3D fields.")

        # Extend the field if not a complete array is given
        field = self.plot_tool_extend_field(field)

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
        self.plot_tool_set_axis_properties(**kwargs)

        return fig, ax

    def plot_complex_field_in_plane(self, complex_field, normal_vector=None, position=None, **kwargs):
        """
        Plots the complex field in a plane perpendicular to the given normal vector using

        Args:
            complex_field (array-like): The complex field to be plotted.
            normal_vector (array-like, optional): The normal vector of the plane. Default is [0,1,0].
            position (array-like, optional): The position of the plane. Default is the middle of the system.
            **kwargs: Keyword arguments for the plot, see https://vidarsko.github.io/ComFiT/Plotting/.

        Returns:
            matplotlib.axes.Axes: The axes containing the plot of the complex field.
        """

        if self.dim != 3:
            raise Exception("The plot in plane function is only defined for 3D fields.")

        # Extend the field if not a complete array is given
        complex_field = self.plot_tool_extend_field(complex_field)

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
        self.plot_tool_set_axis_properties(**kwargs)

        return fig, ax

    def plot_angle_field_in_plane(self, angle_field, normal_vector=None, position=None,**kwargs):
        """
        Plots the angle field in a plane.

        Args:
            angle_field (numpy.ndarray): The angle field to be plotted.
            normal_vector (array-like, optional): The normal vector of the plane. Default is [0,1,0].
            position (array-like, optional): The position of the plane. Default is the middle of the system.
            **kwargs: Keyword arguments for the plot, see https://vidarsko.github.io/ComFiT/Plotting/.

        Returns:
            matplotlib.axes.Axes: The axes containing the plot.
        """

        # Check if the vector field is complex
        if np.iscomplexobj(angle_field):
            print("\033[91mWarning: the angle vector field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
            print('Max imaginary part: ', np.max(np.imag(angle_field)))
            angle_field = np.real(angle_field)

        # Extend the field if not a complete array is given
        angle_field = self.plot_tool_extend_field(angle_field)

        complex_field = np.exp(1j * angle_field)
    
        return self.plot_complex_field_in_plane(complex_field, normal_vector=normal_vector, position=position, **kwargs)

    
    def plot_vector_field_in_plane(self,vector_field,position=None,normal_vector=None,spacing=2,**kwargs):
        """
        Plots the vector field in a plane.
        
        Args:
            vector_field (tuple): Tuple containing the x and y components of the vector field.
            position (array-like, optional): The position of the plane. Default is the middle of the system.
            normal_vector (array-like, optional): The normal vector of the plane. Default is [0,1,0].
            spacing (int, optional): The spacing for the quiver plot. Default is 5.
            kwargs: Keyword arguments for the plot, see https://vidarsko.github.io/ComFiT/Plotting/.

        Returns:
            matplotlib.figure.Figure: The figure containing the plot.
            matplotlib.axes.Axes: The axes containing the plot.
        """

        # Convert the vector field to a numpy array
        vector_field = np.array(vector_field)

        # Extend the field if not a complete array is given
        # TODO: Make this part more efficient (Vidar 11.03.24)
        vector_field_copy = []
        for n in range(vector_field.shape[0]):
            vector_field_copy.append(self.plot_tool_extend_field(vector_field[n]))
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
        self.plot_tool_set_axis_properties(**kwargs)
        return fig, ax