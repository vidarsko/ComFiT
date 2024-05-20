import numpy as np
import scipy as sp

from comfit.tools.tool_colormaps import tool_colormap_angle, tool_colormap_bluewhitered, tool_colormap_sunburst
from comfit.tools.tool_create_orthonormal_triad import tool_create_orthonormal_triad

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.tri as mtri
import matplotlib.cm as cm

import plotly.graph_objects as go

from skimage.measure import marching_cubes

class BaseSystemPlot:
    """ Plotting methods for the base system class"""

    def show(self):
        """Show the plot."""
        if self.plot_lib == 'matplotlib':
            plt.show()
        elif self.plot_lib == 'plotly':
            go.show() #or fig.show?

    def plot_tool_set_axis_properties(self, **kwargs):
        """Sets the properties of the axis for a plot.
        
        Args:
            kwargs: keyword arguments for the axis properties

        Returns:
            The axis object with the properties set.
        """

        plot_lib = kwargs.get('plot_lib', self.plot_lib)
        
        if plot_lib == 'matplotlib':
            ax = kwargs.get('ax', plt.gca())
        elif plot_lib == 'plotly':
            fig = kwargs.get('fig', go.Figure())

        xlabel = kwargs.get('xlabel', 'x/a₀')

        if 'title' in kwargs:
                ax.set_title(kwargs['title'])

        #TODO: Check if I should instead include a fig object somewhere (Vidar 05.02.24)
        suptitle = kwargs.get('suptitle', None)
        if suptitle is not None:
            if plot_lib == 'matplotlib':
                ax.get_figure().suptitle(suptitle)
            elif plot_lib == 'plotly':
                fig.update_layout(title_text=suptitle)
        
        grid = kwargs.get('grid', True)
        if plot_lib == 'matplotlib':
            ax.grid(grid)

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
        xticks = kwargs.get('xticks', None)
        if xticks is not None:
            if plot_lib == 'matplotlib':
                ax.set_xticks(xticks)
        
        xticklabels = kwargs.get('xticklabels', None)
        if xticklabels is not None:
            if plot_lib == 'matplotlib':
                ax.set_xticklabels(xticklabels)
            elif plot_lib == 'plotly':
                fig.update_layout(xaxis=dict(tickvals=xticks, ticktext=xticklabels))

        yticks = kwargs.get('yticks', None)
        if yticks is not None:
            if plot_lib == 'matplotlib':
                ax.set_yticks(yticks)
            elif plot_lib == 'plotly':
                fig.update_layout(yaxis=dict(tickvals=yticks))
        
        yticklabels = kwargs.get('yticklabels', None)
        if yticklabels is not None:
            if plot_lib == 'matplotlib':
                ax.set_yticklabels(yticklabels)
            elif plot_lib == 'plotly':
                fig.update_layout(yaxis=dict(ticktext=yticklabels))

        zticks = kwargs.get('zticks', None)
        if zticks is not None:
            if plot_lib == 'matplotlib':
                ax.set_zticks(zticks)
            elif plot_lib == 'plotly':
                fig.update_layout(zaxis=dict(tickvals=zticks))
        
        zticklabels = kwargs.get('zticklabels', None)
        if zticklabels is not None:
            if plot_lib == 'matplotlib':
                ax.set_zticklabels(zticklabels)
            elif plot_lib == 'plotly':
                fig.update_layout(zaxis=dict(ticktext=zticklabels))
        
        if self.dim == 1:
            # Set the xlabel
            if plot_lib == 'matplotlib':
                ax.set_xlabel(xlabel)
            elif plot_lib == 'plotly':
                fig.update_layout(xaxis_title=xlabel)

            # Set the ylabel
            if 'ylabel' in kwargs:
                if plot_lib == 'matplotlib':
                    ax.set_ylabel(kwargs['ylabel'])
                elif plot_lib == 'plotly':
                    fig.update_layout(yaxis_title=kwargs['ylabel'])

            # Set the title
            if 'title' in kwargs:
                if plot_lib == 'matplotlib':
                    ax.set_title(kwargs['title'])
                elif plot_lib == 'plotly':
                    fig.update_layout(title_text=kwargs['title'])
            
            if plot_lib == 'matplotlib':
                ax.set_xlim(xlim[0], xlim[1])
            elif plot_lib == 'plotly':
                fig.update_layout(xaxis_range=[xlim[0], xlim[1]])

            if ylim is not None:
                if plot_lib == 'matplotlib':
                    ax.set_ylim(ylim[0], ylim[1])
                elif plot_lib == 'plotly':
                    fig.update_layout(yaxis_range=[ylim[0], ylim[1]])

            if zlim is not None:
                if plot_lib == 'matplotlib':
                    ax.set_zlim(zlim[0], zlim[1])
                elif plot_lib == 'plotly':
                    fig.update_layout(zaxis_range=[zlim[0], zlim[1]])
        
        elif self.dim == 2:
            
            ylabel = kwargs.get('ylabel', 'y/a₀')

            if plot_lib == 'matplotlib':
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
            elif plot_lib == 'plotly':
                fig.update_layout(xaxis_title=xlabel)
                fig.update_layout(yaxis_title=ylabel)
            
            if plot_lib == 'matplotlib':
                ax.set_xlim(xlim[0], xlim[1])
                ax.set_ylim(ylim[0], ylim[1])
            elif plot_lib == 'plotly':
                fig.update_layout(xaxis_range=[xlim[0], xlim[1]])
                fig.update_layout(yaxis_range=[ylim[0], ylim[1]])

            if zlim is not None:
                if plot_lib == 'matplotlib':
                    ax.set_zlim(zlim[0], zlim[1])
                elif plot_lib == 'plotly':
                    fig.update_layout(zaxis_range=[zlim[0], zlim[1]])

            # Set the aspect ratio
            axis_equal = kwargs.get('axis_equal', True)
            if axis_equal:
                if plot_lib == 'matplotlib':
                    ax.set_aspect('equal')

        elif self.dim == 3:

            ylabel = kwargs.get('ylabel', 'y/a₀')
            zlabel = kwargs.get('zlabel', 'z/a₀')

            if plot_lib == 'matplotlib':
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.set_zlabel(zlabel)

                ax.set_xlim3d(xlim[0], xlim[1])
                ax.set_ylim3d(ylim[0], ylim[1])
                ax.set_zlim3d(zlim[0], zlim[1])
            elif plot_lib == 'plotly':
                fig.update_layout(
                scene = dict(
                    xaxis=dict(title=xlabel, range=[xlim[0], xlim[1]]),
                    yaxis=dict(title=ylabel, range=[ylim[0], ylim[1]]),
                    zaxis=dict(title=zlabel, range=[zlim[0], zlim[1]])
                    )       
                    )

            elif plot_lib == 'plotly':
                fig.update_layout(xaxis_range=[xlim[0], xlim[1]])
                fig.update_layout(yaxis_range=[ylim[0], ylim[1]])
                fig.update_layout(zaxis_range=[zlim[0], zlim[1]])

            # Set the aspect ratio
            axis_equal = kwargs.get('axis_equal', True)
            if axis_equal:
                if plot_lib == 'matplotlib':
                    ax.set_aspect('equal')  
        if plot_lib == 'matplotlib':
            return ax
        elif plot_lib == 'plotly':
            return fig

    def plot_tool_surface_matplotlib(self, **kwargs):
        """Plots the surface of the given field.

        Args:
            **kwargs: Keyword arguments for the plot.
        
        Returns:
            The axes containing the plot.
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

        return ax

    def plot_tool_extend_field(self,field):
        """Extends the field in case not a complete array is given. 

        For instance, if a field in 3 dimensions is calculated using only x, then the field is extended to the full 3D array.

        Args:
            field: The field to be extended

        Returns:
            The extended field (np.ndarray)
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

    def plot_tool_plotly_add_angle_colorbar3D(self,fig):
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
        """Plots the given (real) field.
        
        Args:
            field (array-like): The field to be plotted.
            **kwargs: Keyword arguments for the plot.
                See https://comfitlib.com/ClassBaseSystem/ 
                for a full list of keyword arguments.
        
        Returns:
            If self.plot_lib == 'matplotlib':
                The axes containing the plot (matplotlib.axes.Axes).
            If self.plot_lib == 'plotly':
                The figure containing the plot.
        """

        if field.dtype == bool:
            field = field.astype(float)

        # Check if the vector field is complex
        if np.iscomplexobj(field):
            print("\033[91mWarning: the provided field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
            print('Max imaginary part: ', np.max(np.imag(field)))
            field = np.real(field)

        plot_lib = kwargs.get('plot_lib', self.plot_lib)

        if plot_lib=='matplotlib':
            # Check if an axis object is provided
            fig = kwargs.get('fig', plt.gcf())
            ax = kwargs.get('ax', None)
        elif plot_lib=='plotly':
            fig = kwargs.get('fig', go.Figure())

        # Kewyord arguments
        colorbar = kwargs.get('colorbar', True)
        axis_equal = kwargs.get('axis_equal',True)

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

            if self.plot_lib == 'matplotlib':
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
            elif self.plot_lib == 'plotly':
                X, Y = np.meshgrid(self.x, self.y, indexing='ij')

                fig.add_trace(go.Heatmap(
                    x=X.flatten()/self.a0,
                    y=Y.flatten()/self.a0,
                    z=field.flatten(),
                    zmin=np.min(field),
                    zmax=np.max(field),
                    colorscale='Viridis'
                ))

                if axis_equal:
                    # Set axis to be equal
                    fig.update_yaxes(
                        scaleanchor="x",
                        scaleratio=1,
                    )

                    fig.update_xaxes(
                        scaleanchor="y",
                        scaleratio=1,
                        ) 

        elif self.dim == 3:

            # Keyword arguments particular to the 3D case

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
            

            if 'colorbar' in kwargs:
                colorbar = kwargs['colorbar']
            else:
                colorbar = True

            #Plotting the layers
            if plot_lib == 'matplotlib':
                if ax == None:
                    plt.clf()
                    ax = plt.gcf().add_subplot(111, projection='3d')


                if field_min < layer_values[1] < field_max:
                    self.plot_tool_surface_matplotlib(field=field, 
                                            value=layer_values[1], 
                                            color=colormap((layer_values[1]-vmin) / (vmax-vmin)), 
                                            alpha=alpha,
                                            ax=ax)

                for layer_value in layer_values[2:-1]:
                    if field_min < layer_value < field_max:
                        self.plot_tool_surface_matplotlib(field=field, 
                                                value=layer_value, 
                                                color=colormap((layer_value-vmin) / (vmax-vmin)), 
                                                alpha=alpha,
                                                ax=ax)


                if colorbar:
                    sm = plt.cm.ScalarMappable(cmap=colormap)
                    sm.set_clim(vmin, vmax)
                    plt.colorbar(sm, ax=ax, pad=0.2)

            elif plot_lib == 'plotly':

                X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

                for layer_value in layer_values[1:-1]:

                    fig.add_trace(go.Isosurface(
                        x=X.flatten()/self.a0,
                        y=Y.flatten()/self.a0,
                        z=Z.flatten()/self.a0,
                        value = field.flatten(),
                        isomin = layer_value,
                        isomax = layer_value,
                        cmin = vmin,
                        cmax = vmax,
                        surface=dict(count=3),  # Ensuring only one surface is shown
                        colorscale='Viridis',
                        opacity=alpha,
                        showscale=bool(layer_value == layer_values[1])
                    )
                        )
        
        if plot_lib == 'matplotlib':
            kwargs['ax'] = ax
            self.plot_tool_set_axis_properties(**kwargs)
            return fig, ax

        elif plot_lib == 'plotly':
            kwargs['fig'] = fig
            self.plot_tool_set_axis_properties(**kwargs)
            return fig

        

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

        plot_lib = kwargs.get('plot_lib', self.plot_lib)

        if plot_lib == 'matplotlib':
            # Check if an axis object is provided
            fig = kwargs.get('fig', plt.gcf())
            ax = kwargs.get('ax', None)
        elif plot_lib == 'plotly':
            fig = kwargs.get('fig', go.Figure())

        # Kewyord arguments
        colorbar = kwargs.get('colorbar', True)
        axis_equal = kwargs.get('axis_equal',True)

        # Extend the field if not a complete array is given
        complex_field = self.plot_tool_extend_field(complex_field)

        # Calculate the magnitude and phase of the complex field
        rho = np.abs(complex_field)
        theta = np.angle(complex_field)

        if self.dim == 1:

            # Keyword arguments particular to the 1D case
            grid = kwargs.get('grid', False)

            if self.plot_lib == 'matplotlib':
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

            elif self.plot_lib == 'plotly':
                
                if fig == None:
                    fig = go.Figure()

                # Color in the graph based on the argument of the complex field
                blend_factor=0.3 # The degree to which the color is blended with white
                cmap = tool_colormap_angle()

                color = (1-blend_factor)*np.array(cmap((theta[0] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1])
                color_str = ('rgb('+str(int(color[0]*255))+','+str(int(color[1]*255))+','+str(int(color[2]*255))+')')

                fig.add_trace(go.Scatter(x=[self.xmin/self.a0,(self.xmin+self.dx/2)/self.a0], 
                              y=[rho[0],(rho[0]+rho[1])/2],
                                mode='lines',
                                line=dict(color='rgba(0,0,0,0)'),
                                fill='tozeroy',
                                showlegend=False,
                                hoverinfo='skip',
                                fillcolor=color_str))

                for i in range(1,self.xRes-1):
                    color = (1-blend_factor)*np.array(cmap((theta[i] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1])
                    color_str = ('rgb('+str(int(color[0]*255))+','+str(int(color[1]*255))+','+str(int(color[2]*255))+')')
                    fig.add_trace(go.Scatter(x=[(self.x[i]-self.dx/2)/self.a0,self.x[i]/self.a0], 
                              y=[(rho[i]+rho[i-1])/2,rho[i]],
                                mode='lines',
                                line=dict(color='rgba(0,0,0,0)'),
                                fill='tozeroy',
                                showlegend=False,
                                hoverinfo='skip',
                                fillcolor=color_str))

                    color = (1-blend_factor)*np.array(cmap((theta[i] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1])
                    color_str = ('rgb('+str(int(color[0]*255))+','+str(int(color[1]*255))+','+str(int(color[2]*255))+')')
                    fig.add_trace(go.Scatter(x=[self.x[i]/self.a0,(self.x[i]+self.dx/2)/self.a0], 
                              y=[rho[i],(rho[i]+rho[i+1])/2],
                                mode='lines',
                                line=dict(color='rgba(0,0,0,0)'),
                                fill='tozeroy',
                                showlegend=False,
                                hoverinfo='skip',
                                fillcolor=color_str))

                color = (1-blend_factor)*np.array(cmap((theta[-1] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1])
                color_str = ('rgb('+str(int(color[0]*255))+','+str(int(color[1]*255))+','+str(int(color[2]*255))+')')
                fig.add_trace(go.Scatter(x=[(self.xmax-1.5*self.dx)/self.a0,(self.xmax-self.dx)/self.a0], 
                              y=[(rho[-1]+rho[-2])/2,rho[-1]],
                                mode='lines',
                                line=dict(color='rgba(0,0,0,0)'),
                                fill='tozeroy',
                                showlegend=False,
                                hoverinfo='skip',
                                fillcolor=color_str))

                fig.add_trace(go.Scatter(
                    x=self.x/self.a0,
                    y=rho,
                    mode='lines',
                    showlegend=False,
                    customdata=np.stack((theta/np.pi, rho), axis=-1),
                    hovertemplate='x: %{x:.2f} a₀<br>θ: %{customdata[0]:.2f} π<br>ρ: %{customdata[1]:.2e}',
                    name='',
                    line=dict(color='black')
                ))
                

        elif self.dim == 2:
            # Keyword arguments particular to the 2D case
            plot_method = kwargs.get('plot_method', 'phase_angle')

            # Create a meshgrid
            X, Y = np.meshgrid(self.x, self.y, indexing='ij')

            if plot_method == '3Dsurface':
                if self.plot_lib == 'matplotlib':  
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

                elif self.plot_lib == 'plotly':
                    print('Plotly not yet implemented for 2D 3Dsurface angle plots.')
                    pass

            elif plot_method == 'phase_angle':
                

                # Keyword arguments particular to the phase angle plot
                grid = kwargs.get('grid', False)

                rho_normalized = rho / np.max(rho)

                custom_colormap = tool_colormap_angle()
                
                if self.plot_lib == 'matplotlib':
                    
                    # Check if an axis object is provided
                    if ax == None:
                        fig.clf()
                        ax = fig.add_subplot(111)

                    # Padding for the colorbar
                    padding=0.05
                    
                    mesh = ax.pcolormesh(X/self.a0, Y/self.a0, theta, shading='auto', cmap=custom_colormap, vmin=-np.pi, vmax=np.pi)
                    mesh.set_alpha(rho_normalized)
                
                elif self.plot_lib == 'plotly':

                    norm = mcolors.Normalize(vmin=-np.pi, vmax=np.pi)
                    colormap = cm.ScalarMappable(norm=norm, cmap=custom_colormap)

                    # Create an RGBA image array
                    rgba_image = colormap.to_rgba(theta)

                    # Apply the spatial opacity
                    for i in range(rgba_image.shape[0]):
                        for j in range(rgba_image.shape[1]):
                            for k in range(3):
                                rgba_image[i, j, k] = 1-(1-rgba_image[i, j, k])*rho_normalized[i, j]  # Set alpha channel

                    rgba_image = np.transpose(rgba_image, (1, 0, 2))

                    # Convert RGBA image to a format Plotly can understand
                    image_data = (rgba_image * 255).astype(np.uint8)


                    if fig == None:
                        fig = go.Figure()

                    fig.add_trace(go.Image(z=image_data, 
                                    dx=self.dx/self.a0, 
                                    dy=self.dy/self.a0, 
                                    x0=self.xmin/self.a0, 
                                    y0=self.ymin/self.a0,
                                    hovertemplate='x: %{x:.2f} a₀<br>y: %{y:.2f} a₀<br>θ: %{customdata[0]:.2f} π<br>ρ: %{customdata[1]:.2e}',
                                    customdata=np.stack((np.transpose(theta/np.pi), np.transpose(rho)), axis=-1),
                                    name=''
                                    ))  


                    # fig.add_trace(go.Heatmap(
                    #     x=X.flatten()/self.a0,
                    #     y=Y.flatten()/self.a0,
                    #     z=theta.flatten(),
                    #     zmin=-np.pi,
                    #     zmax=np.pi,
                    #     colorscale=tool_colormap_angle(pyplot=True),
                    #     colorbar=dict(
                    #             tickvals=[-np.pi, -2*np.pi/3, -np.pi/3, 0, np.pi/3, 2*np.pi/3, np.pi],
                    #             ticktext=['-π', '-2π/3', '-π/3', '0', 'π/3', '2π/3', 'π']
                    #         )
                    # ))

                    if axis_equal:
                       # Set axis to be equal
                        fig.update_yaxes(
                            scaleanchor="x",
                            scaleratio=1,
                        )

                        fig.update_xaxes(
                            scaleanchor="y",
                            scaleratio=1,
                            ) 

                    

                    

        elif self.dim == 3:

            grid = kwargs.get('grid', True)

            plot_method = kwargs.get('plot_method', 'phase_blob')
            
            X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

            rho_normalized = rho / np.max(rho)

            colormap = tool_colormap_angle()

            if plot_lib == 'matplotlib':
                ax = kwargs.get('ax', None)

                if ax == None:
                    fig.clf()
                    ax = fig.add_subplot(111, projection='3d')
            
            elif plot_lib == 'plotly':
                pass
        
            if plot_method == 'phase_angle':
                
                if self.plot_lib == 'matplotlib':
                    
                    # Padding for the colorbar
                    padding=0.2

                    for angle in [-2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3]:
                        field_to_plot = theta.copy()
                        field_to_plot[theta < angle - 1] = float('nan')
                        field_to_plot[theta > angle + 1] = float('nan')
                        field_to_plot[rho_normalized < 0.01] = float('nan')

                        if np.nanmin(field_to_plot) < angle < np.nanmax(field_to_plot):
                            #TODO: make alpha a keyword argument (Vidar 11.03.24)
                            self.plot_tool_surface_matplotlib(field=field_to_plot, 
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
                        self.plot_tool_surface_matplotlib(field=field_to_plot, 
                                                value=np.pi, 
                                                color=colormap(0), 
                                                alpha=0.5,
                                                ax=ax)
                elif self.plot_lib == 'plotly':
                    print('Plotly not yet implemented for 3D phase angle plots.')
                    pass
            
            elif plot_method == 'phase_blob':
                # TODO: change the function so that it used plot_tool_surface_matplotlib (Vidar 11.03.24)
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

                    if plot_lib == 'matplotlib':
                        

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

                    elif plot_lib == 'plotly':
                        # Convert colors to 'rgba()' format required by Plotly
                        colors = ['rgba({},{},{},{})'.format(*c) for c in (255*colors).astype(int)]

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
                            facecolor=colors,  # Set color for each face
                            showscale=True
                        )

                        fig = go.Figure(data=[mesh])

        if colorbar:
            if plot_lib == 'matplotlib':
                mappable = plt.cm.ScalarMappable(cmap=tool_colormap_angle())
                mappable.set_array([])
                mappable.set_clim(-np.pi, np.pi)
                cbar = plt.colorbar(mappable, ax=ax, pad=padding)
                cbar.set_ticks(np.array([-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi]))
                cbar.set_ticklabels([r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$'])

            elif plot_lib == 'plotly':
                # Add a colorbar to the plot
                if self.dim == 1:
                    pass
                elif self.dim == 2:
                    pass
                elif self.dim == 3:
                    fig = self.plot_tool_plotly_add_angle_colorbar3D(fig)

        if plot_lib == 'matplotlib':
            kwargs['ax'] = ax
            kwargs['grid'] = grid
            self.plot_tool_set_axis_properties(**kwargs)
            return fig, ax

        elif plot_lib == 'plotly':
            kwargs['fig'] = fig
            self.plot_tool_set_axis_properties(**kwargs)
            return fig

    def plot_angle_field(self, angle_field, **kwargs):
        """Plot the angle field.

        Args:
            field (array-like): The angle field values.
            ax (matplotlib.axes.Axes, optional): The axes to plot the angle field on. If not provided, a new subplot will be created.
        
        Returns:
            The axes containing the plot. (matplotlib.axes.Axes)
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
        """Plots a vector field.

        Args:
            vector_field (tuple): Tuple containing the x and y components of the vector field.
            spacing (int, optional): The spacing for the quiver plot. Default is 5.
            **kwargs: Keyword arguments for the plot, see https://comfitlib.com/Plotting/.

        Returns:
            Tuple containing
                - The figure containing the plot.
                - The axes containing the plot.
        """
        # Convert the vector field to a numpy array
        vector_field = np.array(vector_field)
        
        # Extend the field if not a complete array is given
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
                    fig.clf()
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
                    fig.clf()
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
                    fig.clf()
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