from typing import TYPE_CHECKING, Any, Dict, Tuple

if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

# Standard library imports
import numpy as np
import scipy as sp

# Third-party library imports
import plotly.graph_objects as go
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from skimage.measure import marching_cubes

from .plot_surface_plotly import plot_surface_plotly

# Local application imports
from comfit.tool import (
    tool_complete_field,
    tool_colormap, 
    tool_set_plot_axis_properties_plotly,
    tool_plotly_colorbar,
    tool_plotly_define_2D_plot_ax,
    tool_plotly_define_3D_plot_ax
)


def plot_complex_field_plotly(
        self: 'BaseSystem',
        complex_field: np.ndarray,
        **kwargs: Any
        ) -> Tuple[go.Figure, Dict]:
    """Plot a complex field using Plotly.

    Creates a visualization of a complex field using various plotting methods
    depending on the dimensionality of the field.

    Parameters
    ----------
    self : BaseSystem
        A BaseSystem (or derived) instance.
    complex_field : np.ndarray
        The complex field to plot
    kwargs : Any
        Additional keyword arguments to customize the plot

    Returns
    -------
    go.Figure
        The Plotly figure containing the plot
    """

    kwargs['colormap'] = kwargs.get('colormap', 'angle') # Override the default colormap with 'angle'
    complex_field, fig, ax, kwargs = self.plot_prepare(complex_field, field_type = 'complex', **kwargs)

    # Kewyord arguments
    kwargs['colorbar'] = kwargs.get('colorbar', True)

    # Extend the field if not a complete array is given
    complex_field = tool_complete_field(self, complex_field)

    # Calculate the magnitude and phase of the complex field
    rho = np.abs(complex_field)
    theta = np.angle(complex_field)

    plt_colormap_object = tool_colormap(kwargs['colormap'], plot_lib='matplotlib')

    ###############################################################
    ###################### DIMENSION: 1 ###########################
    ###############################################################

    if self.dim == 1:

        ax = tool_plotly_define_2D_plot_ax(fig, ax)

        vlim = kwargs.get('vlim', None)
        if vlim is not None:
            kwargs['ylim'] = vlim

        # Color in the graph based on the argument of the complex field
        blend_factor=0.3 # The degree to which the color is blended with white
        
        color = (1-blend_factor)*np.array(plt_colormap_object((theta[0] + np.pi) / (2 * np.pi)))+blend_factor*np.array([1,1,1,1])
        color_str = ('rgb('+str(int(color[0]*255))+','+str(int(color[1]*255))+','+str(int(color[2]*255))+')')

        # Extract coordinates
        x = kwargs.get('x', self.x/self.a0).flatten()
        dx = x[1] - x[0]
        xmin = x[0]
        xmax = x[-1]+dx

        fig.add_trace(go.Scatter(x=[xmin,(xmin+dx/2)], 
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
            fig.add_trace(go.Scatter(x=[(x[i]-dx/2),x[i]], 
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
            fig.add_trace(go.Scatter(x=[x[i],(x[i]+dx/2)], 
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
        fig.add_trace(go.Scatter(x=[(xmax-1.5*dx),(xmax-dx)], 
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
            x=x,
            y=rho,
            mode='lines',
            showlegend=False,
            customdata=np.stack((theta/np.pi, rho), axis=-1),
            hovertemplate= kwargs['xlabel']+': %{x:.2f} <br>'+\
                                'amplitude: %{customdata[1]:.2e}<br>'+\
                                'phase: %{customdata[0]:.2f} π',
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

        # Extract coordinates
        x = kwargs.get('x', self.x/self.a0).flatten()
        dx = x[1] - x[0]
        xmin = x[0]
        xmax = x[-1]+dx
        
        y = kwargs.get('y', self.y/self.a0).flatten()
        dy = y[1] - y[0]
        ymin = y[0]
        ymax = y[-1]+dy

        rho_normalized = rho / np.max(rho)

        if plot_method == '3Dsurface':

            ax = tool_plotly_define_3D_plot_ax(fig, ax)

            # Plot field as 3D surface where height is amplitude and color is phase
            X, Y = np.meshgrid(x, y, indexing='ij')

            surface = go.Surface(
                x=X,
                y=Y,
                z=rho,
                surfacecolor=theta,
                colorscale=kwargs['colormap_object'],
                customdata=np.stack((theta/np.pi, rho), axis=-1),
                hovertemplate=(f"{kwargs['xlabel']}: %{{x:.2f}}<br>" +
                              f"{kwargs['ylabel']}: %{{y:.2f}}<br>" +
                              "amplitude: %{customdata[1]:.2e}<br>" +
                              "phase: %{customdata[0]:.2f} π"),
                showscale=False
            )

            fig.add_trace(surface)


        elif plot_method == 'phase_angle':
            
            ax = tool_plotly_define_2D_plot_ax(fig, ax)
            
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

            trace = go.Image(z=image_data, 
                            dx=dx, 
                            dy=dy, 
                            x0=xmin, 
                            y0=ymin,
                            hovertemplate=kwargs['xlabel']+': %{x:.2f}<br>'+\
                                            kwargs['ylabel']+': %{y:.2f}<br>'+\
                                            'amplitude: %{customdata[1]:.2e}<br>'+\
                                            'phase: %{customdata[0]:.2f} π',
                            customdata=np.stack((np.transpose(theta/np.pi), np.transpose(rho)), axis=-1),
                            name='',
                            xaxis=ax['xN'],
                            yaxis=ax['yN']
                            )

            fig.add_trace(trace)

            # Add grid lines manually since the default plotly grid lines are
            # not available for go.Image traces
            if kwargs['grid']:
                N = 11
                # Draw N equally spaced horizontal and vertical lines
                x_positions = np.linspace(x[0], x[-1], N)
                y_positions = np.linspace(y[0], y[-1], N)

                # Add horizontal grid lines
                for y_pos in y_positions:
                    fig.add_trace(go.Scatter(
                        x=[xmin, xmax],
                        y=[y_pos, y_pos],
                        mode='lines',
                        line=dict(color='grey', width=0.5),
                        showlegend=False,
                        hoverinfo='skip',
                        xaxis=ax['xN'],
                        yaxis=ax['yN']
                    ))

                # Add vertical grid lines
                for x_pos in x_positions:
                    fig.add_trace(go.Scatter(
                        x=[x_pos, x_pos],
                        y=[ymin, ymax],
                        mode='lines',
                        line=dict(color='grey', width=0.5),
                        showlegend=False,
                        hoverinfo='skip',
                        xaxis=ax['xN'],
                        yaxis=ax['yN']
                    ))


    ###############################################################
    ###################### DIMENSION: 3 ###########################
    ###############################################################
    elif self.dim == 3:

        ax = tool_plotly_define_3D_plot_ax(fig, ax)
        kwargs['ax'] = ax 

        # Keyword arguments particular to the 3D case
        plot_method = kwargs.get('plot_method', 'phase_blob')
        
        rho_normalized = rho / np.max(rho)

        if plot_method == 'phase_angle':

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

            X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
            
            for angle in [-2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3]:
                field_to_plot = theta.copy()
                field_to_plot[theta < angle - 1] = float('nan')
                field_to_plot[theta > angle + 1] = float('nan')
                field_to_plot[rho_normalized < 0.01] = float('nan')

                if np.nanmin(field_to_plot) < angle < np.nanmax(field_to_plot):

                    mesh, shadow_meshes = plot_surface_plotly(self,
                                            field_to_plot,
                                            value=angle,
                                            alpha=rho,
                                            color=plt_colormap_object((angle + np.pi) / (2 * np.pi)), 
                                            plt_colormap_object=plt_colormap_object,
                                            **kwargs
                                            )

                    fig.add_trace(mesh)

            theta = np.mod(theta, 2 * np.pi)

            field_to_plot = theta.copy()
            field_to_plot[theta < np.pi - 1] = float('nan')
            field_to_plot[theta > np.pi + 1] = float('nan')
            field_to_plot[rho_normalized < 0.01] = float('nan')

            if np.nanmin(field_to_plot) < np.pi < np.nanmax(field_to_plot):
                mesh, shadow_meshes = plot_surface_plotly(self,
                                            field_to_plot,
                                            value=np.pi,
                                            alpha=rho,
                                            color=plt_colormap_object((np.pi + np.pi) / (2 * np.pi)), 
                                            plt_colormap_object=plt_colormap_object,
                                            **kwargs
                                            )

                fig.add_trace(mesh)
        
        elif plot_method == 'phase_blob':
            
            phase_blob_threshold = kwargs.get('phase_blob_threshold', 0.5)

            if np.nanmin(rho_normalized)<phase_blob_threshold<np.nanmax(rho_normalized):

                mesh, shadow_meshes = plot_surface_plotly(self, 
                                        field=rho, 
                                        value=phase_blob_threshold*np.max(rho), 
                                        alpha=1,
                                        color=complex_field,
                                        plt_colormap_object=plt_colormap_object,
                                          **kwargs)

                fig.add_trace(mesh)

                if shadow_meshes is not None:
                    for shadow_mesh in shadow_meshes:
                        fig.add_trace(shadow_mesh)



    # if kwargs['colorbar'] and not(ax['colorbar']) and not(kwargs['field_is_nan']):
    #     ax['colormap_object'] = kwargs['colormap_object']
    #     fig.add_trace(tool_plotly_colorbar(ax, type='angle'))
    #     ax['colorbar'] = True


    if kwargs['colorbar'] and not(kwargs['field_is_nan']):
        ax['colormap_object'] = kwargs['colormap_object']
        fig.add_trace(tool_plotly_colorbar(ax, type='angle'))
        ax['colorbar'] = True

    kwargs['fig'] = fig
    kwargs['ax'] = ax
    tool_set_plot_axis_properties_plotly(self, **kwargs)

    return fig, ax