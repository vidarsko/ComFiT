import numpy as np
import plotly.graph_objects as go
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import scipy as sp
from skimage.measure import marching_cubes
from comfit.tool import tool_complete_field
from comfit.tool import tool_colormap
from comfit.tool import tool_set_plot_axis_properties_plotly
from comfit.tool import tool_plotly_colorbar

from comfit.tool import tool_plotly_define_2D_plot_ax
from comfit.tool import tool_plotly_define_3D_plot_ax


def plot_complex_field_plotly(self, complex_field: np.ndarray, **kwargs) -> go.Figure:
    """
    Plot a complex field.

    Args:
        complex_field (numpy.ndarray): The complex field to plot.
        ax (matplotlib.axes.Axes, optional): The matplotlib axes on which to plot the field.
            If not provided, a new 3D axes will be created.
    
    Returns:
        go.Figure: The figure containing the plot.
    """

    kwargs['colormap'] = kwargs.get('colormap', 'angle') # Override the default colormap with 'angle'
    complex_field, fig, ax, kwargs = self.plot_prepare(complex_field, field_type = 'complex', **kwargs)

    # Calculate the magnitude and phase of the complex field
    rho = np.abs(complex_field)
    theta = np.angle(complex_field)

    plt_colormap_object = tool_colormap(kwargs['colormap'], plot_lib='matplotlib')

    ###############################################################
    ###################### DIMENSION: 1 ###########################
    ###############################################################

    if self.dim == 1:

        ax = tool_plotly_define_2D_plot_ax(fig, ax)

        # Keyword arguments particular to the 1D case
        grid = kwargs.get('grid', False)

        vlim = kwargs.get('vlim', None)
        if vlim is not None:
            kwargs['ylim'] = vlim

        # Color in the graph based on the argument of the complex field
        blend_factor=0.3 # The degree to which the color is blended with white
        
        color = (1-blend_factor)*np.array(plt_colormap_object((theta[0] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1])
        color_str = ('rgb('+str(int(color[0]*255))+','+str(int(color[1]*255))+','+str(int(color[2]*255))+')')

        fig.add_trace(go.Scatter(x=[self.xmin/self.a0,(self.xmin+self.dx/2)/self.a0], 
                        y=[rho[0],(rho[0]+rho[1])/2],
                        mode='lines',
                        line=dict(color='rgba(0,0,0,0)'),
                        fill='tozeroy',
                        showlegend=False,
                        hoverinfo='skip',
                        fillcolor=color_str,
                        xaxis=ax['xN'],
                        yaxis=ax['yN']))

        for i in range(1,self.xRes-1):
            color = (1-blend_factor)*np.array(plt_colormap_object((theta[i] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1])
            color_str = ('rgb('+str(int(color[0]*255))+','+str(int(color[1]*255))+','+str(int(color[2]*255))+')')
            fig.add_trace(go.Scatter(x=[(self.x[i]-self.dx/2)/self.a0,self.x[i]/self.a0], 
                        y=[(rho[i]+rho[i-1])/2,rho[i]],
                        mode='lines',
                        line=dict(color='rgba(0,0,0,0)'),
                        fill='tozeroy',
                        showlegend=False,
                        hoverinfo='skip',
                        fillcolor=color_str, 
                        xaxis=ax['xN'],
                        yaxis=ax['yN']))

            color = (1-blend_factor)*np.array(plt_colormap_object((theta[i] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1])
            color_str = ('rgb('+str(int(color[0]*255))+','+str(int(color[1]*255))+','+str(int(color[2]*255))+')')
            fig.add_trace(go.Scatter(x=[self.x[i]/self.a0,(self.x[i]+self.dx/2)/self.a0], 
                        y=[rho[i],(rho[i]+rho[i+1])/2],
                        mode='lines',
                        line=dict(color='rgba(0,0,0,0)'),
                        fill='tozeroy',
                        showlegend=False,
                        hoverinfo='skip',
                        fillcolor=color_str,
                        xaxis=ax['xN'],
                        yaxis=ax['yN']))

        color = (1-blend_factor)*np.array(plt_colormap_object((theta[-1] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1])
        color_str = ('rgb('+str(int(color[0]*255))+','+str(int(color[1]*255))+','+str(int(color[2]*255))+')')
        fig.add_trace(go.Scatter(x=[(self.xmax-1.5*self.dx)/self.a0,(self.xmax-self.dx)/self.a0], 
                        y=[(rho[-1]+rho[-2])/2,rho[-1]],
                        mode='lines',
                        line=dict(color='rgba(0,0,0,0)'),
                        fill='tozeroy',
                        showlegend=False,
                        hoverinfo='skip',
                        fillcolor=color_str,
                        xaxis=ax['xN'],
                        yaxis=ax['yN']))


        fig.add_trace(go.Scatter(
            x=self.x/self.a0,
            y=rho,
            mode='lines',
            showlegend=False,
            customdata=np.stack((theta/np.pi, rho), axis=-1),
            hovertemplate='x: %{x:.2f} a₀<br>θ: %{customdata[0]:.2f} π<br>ρ: %{customdata[1]:.2e}',
            name='',
            line=dict(color='black'),
            xaxis=ax['xN'],
            yaxis=ax['yN']
        ))

        
    ###############################################################
    ###################### DIMENSION: 2 ###########################
    ###############################################################

    elif self.dim == 2:
        
        # Keyword arguments particular to the 2D case
        plot_method = kwargs.get('plot_method', 'phase_angle')

        # Create a meshgrid
        X, Y = np.meshgrid(self.x, self.y, indexing='ij')

        if plot_method == '3Dsurface':
        
            print("\033[91mWarning: 3D surface plot not yet implemented for Plotly.\033[0m")
            pass

        elif plot_method == 'phase_angle':
            
            ax = tool_plotly_define_2D_plot_ax(fig, ax)

            # Keyword arguments particular to the phase angle plot
            grid = kwargs.get('grid', False)

            rho_normalized = rho / np.max(rho)
            
            norm = mcolors.Normalize(vmin=-np.pi, vmax=np.pi)

            # Create an RGBA image array
            rgba_image = cm.ScalarMappable(norm=norm, cmap=plt_colormap_object).to_rgba(theta)

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

            trace = go.Image(z=image_data, 
                            dx=self.dx/self.a0, 
                            dy=self.dy/self.a0, 
                            x0=self.xmin/self.a0, 
                            y0=self.ymin/self.a0,
                            hovertemplate='x: %{x:.2f} a₀<br>y: %{y:.2f} a₀<br>θ: %{customdata[0]:.2f} π<br>ρ: %{customdata[1]:.2e}',
                            customdata=np.stack((np.transpose(theta/np.pi), np.transpose(rho)), axis=-1),
                            name='',
                            xaxis=ax['xN'],
                            yaxis=ax['yN']
                            )

            fig.add_trace(trace)


    ###############################################################
    ###################### DIMENSION: 3 ###########################
    ###############################################################
    elif self.dim == 3:

        ax = tool_plotly_define_3D_plot_ax(fig, ax)

        # Keyword arguments particular to the 3D case
        plot_method = kwargs.get('plot_method', 'phase_blob')
        
        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')
        rho_normalized = rho / np.max(rho)

        if plot_method == 'phase_angle':
        
            print('\033[91mWarning:\033[0m Plotly not yet implemented for 3D phase angle plots.')
            pass
        
        elif plot_method == 'phase_blob':
            
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
                reals = np.real(complex_field)
                imags = np.imag(complex_field)

                interpolation_method = kwargs.get('interpolation_method', 'linear')
                print("Interpolating points with method ", interpolation_method, ".")
                print("If this process is slow, consider passing 'interpolation_method='nearest' with the plot_complex_field function.")
                print("That will speed up the process, but the plot may look less smooth.")
                reals_faces = sp.interpolate.griddata(points, reals.ravel(), centroids, method='nearest')
                imags_faces = sp.interpolate.griddata(points, imags.ravel(), centroids, method='nearest')
                print("Interpolation done.")
                
                theta_faces = np.arctan2(imags_faces, reals_faces)

                # Normalize theta values for color mapping
                theta_faces_normalized = (theta_faces + np.pi) / (2*np.pi)

                # Map normalized theta values to colors
                colors = plt_colormap_object(theta_faces_normalized)

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
                    showscale=True,
                    scene=ax['sceneN']
                )

                fig.add_trace(mesh)

    
    if kwargs['colorbar'] and not(ax['colorbar']) and not(kwargs['field_is_nan']):
        ax['colormap_object'] = kwargs['colormap_object']
        fig.add_trace(tool_plotly_colorbar(ax, type='angle'))
        ax['colorbar'] = True

    kwargs['fig'] = fig
    kwargs['ax'] = ax
    tool_set_plot_axis_properties_plotly(self, **kwargs)

    return fig, ax