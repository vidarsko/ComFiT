from typing import TYPE_CHECKING, Any, Tuple

if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

# Standard library imports
from typing import Optional, Union

# Third-party imports
import numpy as np
import plotly.graph_objects as go

# Local imports
from comfit.tool.tool_set_plot_axis_properties_plotly import tool_set_plot_axis_properties_plotly
from comfit.tool import (
    tool_complete_field,
    tool_add_spacing_2D,
    tool_add_spacing_3D,
    tool_plotly_define_2D_plot_ax,
    tool_plotly_define_3D_plot_ax,
    tool_plotly_colorbar
)

def plot_vector_field_plotly(
    self: 'BaseSystem',
    vector_field: np.ndarray,
    **kwargs: Any
) -> Tuple[go.Figure, dict]:
    """Plot a vector field using Plotly.

    Parameters
    ----------
    vector_field : np.ndarray
        Array containing the vector field components. Shape depends on dimensionality.
    \*\*kwargs : Any
        spacing : int, optional
            The spacing for the quiver plot. Defaults to max(self.xRes//20, 1).
        Additional keyword arguments for plot customization, see https://comfitlib.com/Plotting/

    Returns
    -------
    Tuple[go.Figure, dict]
        The figure containing the plot and the axes dictionary containing plot properties.

    Raises
    ------
    Exception
        If an invalid field is provided to the plot_vector_field function.
    """
    
    spacing = kwargs.get('spacing', max(self.xRes//20,1))

    vector_field, fig, ax, kwargs = self.plot_prepare(vector_field, field_type = 'vector', **kwargs)

    ###############################################################
    ########### DIMENSION: 1 - VECTOR-DIMENSION: 1 ################
    ###############################################################
        
    if self.dim == 1 and vector_field.shape == (1,self.xRes):

        kwargs['colorbar'] = False
        ax = tool_plotly_define_2D_plot_ax(fig, ax)
        
        X, Y = np.meshgrid(self.x, np.array([0]), indexing='ij')
        
        U = np.zeros(X.shape)
        V = np.zeros(X.shape)

        V[:,0] = vector_field[0]

        X,Y,U,V = tool_add_spacing_2D(X,Y,U,V,spacing)
        
        fig, ax = self.plot_field(vector_field[0], **kwargs)

        
    ###############################################################
    ########### DIMENSION: 1 - VECTOR-DIMENSION: 2 ################
    ###############################################################

    elif self.dim == 1 and vector_field.shape == (2,self.xRes):

        ax = tool_plotly_define_3D_plot_ax(fig, ax)
            
        X, Y, Z = np.meshgrid(self.x, np.array([0]), np.array([0]), indexing='ij')

        U = np.zeros(X.shape)
        V = np.zeros(X.shape)
        W = np.zeros(X.shape)

        V[:,0,0] = vector_field[0]
        W[:,0,0] = vector_field[1]

        ax['vmin'] = kwargs.get('vmin', 0)
        ax['vmax'] = kwargs.get('vmax', np.max(np.sqrt(vector_field[0]**2 + vector_field[1]**2)))

        X,Y,Z,U,V,W = tool_add_spacing_3D(X,Y,Z,U,V,W,spacing)

        if not kwargs['field_is_nan']:
            fig.add_trace(go.Cone(x=X.flatten()/self.a0, 
                                    y=Y.flatten()/self.a0, 
                                    z=Z.flatten()/self.a0, 
                                    u=U.flatten(), 
                                    v=V.flatten(), 
                                    w=W.flatten(), 
                                    colorscale=kwargs['colormap_object'], 
                                    sizemode='scaled', 
                                    sizeref=1, 
                                    showscale=False,
                                    scene = ax['sceneN'],
                                    hovertemplate='x: %{x:.1f} a₀ <br>' +
                                                'ux: %{u:.1f} <br>' +
                                                'uy: %{v:.1f} <br>' +
                                                'uz: %{w:.1f} <br>' +
                                                '|u|: %{customdata:.1f}<extra></extra>',
                                                customdata=np.sqrt(U**2 + V**2 + W**2).flatten()))



        # kwargs['ylim'] = [np.min(vector_field[0]), np.max(vector_field[0])]
        # delta_y = kwargs['ylim'][1] - kwargs['ylim'][0]

        # kwargs['zlim'] = [np.min(vector_field[1]), np.max(vector_field[1])]
        # delta_z = kwargs['zlim'][1] - kwargs['zlim'][0]

        # if delta_y < 0.15*delta_z:
        #     kwargs['ylim'] = [kwargs['ylim'][0] - 0.15*delta_z, kwargs['ylim'][1] + 0.15*delta_z]
        
        # if delta_z < 0.15*delta_y:
        #     kwargs['zlim'] = [kwargs['zlim'][0] - 0.15*delta_y, kwargs['zlim'][1] + 0.15*delta_y]

    ###############################################################
    ########### DIMENSION: 1 - VECTOR-DIMENSION: 3 ################
    ###############################################################

    elif self.dim == 1 and vector_field.shape == (3,self.xRes):
        

        ax = tool_plotly_define_3D_plot_ax(fig, ax)
        
        X, Y, Z = np.meshgrid(self.x, np.array([0]), np.array([0]), indexing='ij')

        U = np.zeros(X.shape)
        V = np.zeros(X.shape)
        W = np.zeros(X.shape)

        # Add the vector field
        U[:,0,0] = vector_field[0]
        V[:,0,0] = vector_field[1]
        W[:,0,0] = vector_field[2]

        # Add spacing
        X,Y,Z,U,V,W = tool_add_spacing_3D(X,Y,Z,U,V,W,spacing)

        # Normalize the vectors
        max_vector = np.max(np.sqrt(U ** 2 + V ** 2 + W ** 2))

        ax['vmin'] = kwargs.get('vmin', 0)
        ax['vmax'] = kwargs.get('vmax', np.max(max_vector))

        # Normalizing
        U = U/max_vector
        V = V/max_vector
        W = W/max_vector

        # Scale factors
        vx_scale = kwargs.get('vx_scale', spacing*self.size_x/15/self.a0)
        vy_scale = kwargs.get('vy_scale', 1)
        vz_scale = kwargs.get('vz_scale', 1)

        # Scaling
        U = vx_scale*U
        V = vy_scale*V
        W = vz_scale*W

        if not kwargs['field_is_nan']:
            fig.add_trace(go.Cone(x=X.flatten()/self.a0, 
                            y=Y.flatten()/self.a0, 
                            z=Z.flatten()/self.a0, 
                            u=U.flatten(), 
                            v=V.flatten(), 
                            w=W.flatten(), 
                            colorscale=kwargs['colormap_object'], 
                            sizemode='scaled', 
                            sizeref=1, 
                            showscale=False, 
                            scene = ax['sceneN'],
                            hovertemplate='<b>x:</b> %{x:.1f} a₀ <br>' +
                                        '<b>ux:</b> %{u:.1f} <br>' +
                                        '<b>uy:</b> %{v:.1f} <br>' +
                                        '<b>uz:</b> %{w:.1f} <br>' +
                                        '<b>|u|:</b> %{customdata:.1f}<extra></extra>',
                                        customdata=np.sqrt(U**2 + V**2 + W**2).flatten()))

        # kwargs['axis_equal'] = False
        # kwargs['ylim'] = [-1,1]
        # kwargs['zlim'] = [-1,1]

    ###############################################################
    ########### DIMENSION: 2 - VECTOR-DIMENSION: 1 ################
    ###############################################################

    elif self.dim == 2 and vector_field.shape == (1,self.xRes,self.yRes):
  
        ax = tool_plotly_define_2D_plot_ax(fig, ax)

        X, Y = np.meshgrid(self.x, self.y, indexing='ij')

        U = vector_field[0]
        V = np.zeros(X.shape)

        X,Y,U,V = tool_add_spacing_2D(X,Y,U,V,spacing)

        # Normalize the vectors
        max_vector = np.max(np.sqrt(U ** 2))
        
        ax['vmin'] = kwargs.get('vmin', 0)
        ax['vmax'] = kwargs.get('vmax', np.max(max_vector))

        U = U / max_vector

        # Scale factors
        vx_scale = kwargs.get('vx_scale', spacing*self.dx/self.a0)

        # Scaling
        U = vx_scale*U

        
        fig = kwargs.get('fig', go.Figure())

        u = U.flatten()
        v = V.flatten()
        
        magnitude = np.sqrt(u**2 + v**2)
        magnitude_normalized = magnitude/max(np.max(magnitude),1e-12)

        angle = np.arctan2(v, u)
        direction = np.array([np.cos(angle), np.sin(angle)]).T
        
        fig.add_trace(go.Scatter(
        x=X.flatten()/self.a0,
        y=Y.flatten()/self.a0,
        xaxis=ax['xN'],
        yaxis=ax['yN'],
        mode='markers',
        showlegend=False,
        marker=dict(symbol='arrow', 
            angle=90-angle.flatten()*180/np.pi, 
            size=2*spacing*magnitude_normalized.flatten(), 
            sizemode='diameter',
            color=magnitude.flatten(), 
            colorscale=kwargs['colormap_object'], 
            showscale=False,
            line=dict(color='black')
            ),
            hovertemplate='<b>x:</b> %{x:.1f}a0<br>' +
                        '<b>y:</b> %{y:.1f}a0<br>' +
                        '<b>ux:</b> %{customdata[0]:.2f}<br>' +
                        '<b>uy:</b> %{customdata[1]:.2f}<br>' +
                        '<b>|u|:</b> %{customdata[2]:.2f}<extra></extra>',
            customdata=np.stack((u.flatten(), v.flatten(),magnitude.flatten() ), axis=-1)  # Adding ux, uy and u as customdata
        )
        )

    ###############################################################
    ########### DIMENSION: 2 - VECTOR-DIMENSION: 2 ################
    ###############################################################

    elif self.dim == 2 and vector_field.shape == (2,self.xRes,self.yRes):

        ax = tool_plotly_define_2D_plot_ax(fig, ax)

        X, Y = np.meshgrid(self.x, self.y, indexing='ij')

        X, Y, U, V = tool_add_spacing_2D(X,Y,vector_field[0],vector_field[1],spacing)

        u = U.flatten()
        v = V.flatten()
        
        magnitude = np.sqrt(u**2 + v**2)

        ax['vmin'] = kwargs.get('vmin', 0)
        ax['vmax'] = kwargs.get('vmax', np.max(magnitude))

        magnitude_max = max(np.max(magnitude),1e-12)
        magnitude_normalized = magnitude/magnitude_max

        angle = np.arctan2(v, u)
        direction = np.array([np.cos(angle), np.sin(angle)]).T
        
        if not kwargs['field_is_nan']:
            fig.add_trace(go.Scatter(
            x=X.flatten()/self.a0,
            y=Y.flatten()/self.a0,
            mode='markers',
            xaxis=ax['xN'],
            yaxis=ax['yN'],
            showlegend=False,
            marker=dict(symbol='arrow', 
                angle=90-angle.flatten()*180/np.pi, 
                size=4*spacing*magnitude_normalized.flatten(), 
                sizemode='diameter',
                color=magnitude.flatten(), 
                colorscale=kwargs['colormap_object'], 
                showscale=False,
                cmin=0,
                cmax=magnitude_max,
                line=dict(color='black')
                ),
                hovertemplate='<b>x:</b> %{x:.2f}a0<br>' +
                    '<b>y:</b> %{y:.2f}a0<br>' +
                    '<b>ux:</b> %{customdata[0]:.2e}<br>' +  
                    '<b>uy:</b> %{customdata[1]:.2e}<br>' +
                    '<b>|u|:</b> %{customdata[2]:.2e}<extra></extra>',
                customdata=np.stack((u.flatten(), v.flatten(), magnitude.flatten()), axis=-1)  # Adding ux, uy and u as customdata
            )
            )


    ###############################################################
    ########### DIMENSION: 2 - VECTOR-DIMENSION: 3 ################
    ###############################################################

    elif self.dim == 2 and vector_field.shape == (3,self.xRes,self.yRes):

        ax = tool_plotly_define_3D_plot_ax(fig, ax)
        
        X, Y, Z = np.meshgrid(self.x, self.y, np.array([0]), indexing='ij')
        U = np.zeros(X.shape)
        V = np.zeros(X.shape)
        W = np.zeros(X.shape)
        
        U[:,:,0] = vector_field[0]
        V[:,:,0] = vector_field[1]
        W[:,:,0] = vector_field[2]

        X,Y,Z,U,V,W = tool_add_spacing_3D(X,Y,Z,U,V,W,spacing)

        max_vector = np.max(np.sqrt(U ** 2 + V ** 2 + W ** 2))

        ax['vmin'] = kwargs.get('vmin', 0)
        ax['vmax'] = kwargs.get('vmax', np.max(max_vector))

        # Normalizing
        U = U / max_vector
        V = V / max_vector
        W = W / max_vector
        
        # Scale factors
        vx_scale = kwargs.get('vx_scale', 2*spacing*self.size_x/max_vector/self.a0)
        vy_scale = kwargs.get('vy_scale', 2*spacing*self.size_y/max_vector/self.a0)
        vz_scale = kwargs.get('vz_scale', spacing/self.a0)

        # Scaling
        U = vx_scale*U
        V = vy_scale*V
        W = vz_scale*W

        if not kwargs['field_is_nan']:
            fig.add_trace(go.Cone(x=X.flatten()/self.a0, 
                                y=Y.flatten()/self.a0, 
                                z=Z.flatten()/self.a0, 
                                u=U.flatten(), 
                                v=V.flatten(), 
                                w=W.flatten(), 
                                colorscale=kwargs['colormap_object'], 
                                sizemode='scaled', 
                                sizeref=1, 
                                showscale=False,
                                scene=ax['sceneN'],
                                hovertemplate='<b>x:</b> %{x:.1f} a₀ <br>' +
                                                '<b>y:</b> %{y:.1f} a₀ <br>' +
                                                '<b>z:</b> %{z:.1f} a₀ <br>' +
                                                '<b>ux:</b> %{u:.1f} <br>' +
                                                '<b>uy:</b> %{v:.1f} <br>' +
                                                '<b>uz:</b> %{w:.1f} <br>' +
                                                '<b>|u|:</b> %{customdata:.1f}<extra></extra>',
                                                customdata=np.sqrt(U**2 + V**2 + W**2).flatten()))



        kwargs['axis_equal'] = False
        kwargs['zlim'] = [-spacing,spacing]


    ###############################################################
    ########### DIMENSION: 3 - VECTOR-DIMENSION: 1 ################
    ###############################################################

    elif self.dim == 3 and vector_field.shape == (1,self.xRes,self.yRes,self.zRes):

        ax = tool_plotly_define_3D_plot_ax(fig, ax)

        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')              

        # Define the vector field
        U = vector_field[0]
        V = np.zeros(U.shape)
        W = np.zeros(U.shape)

        # Add spacing
        X,Y,Z,U,V,W = tool_add_spacing_3D(X,Y,Z,U,V,W,spacing)

        # Normalize the vectors
        max_vector = np.max(np.sqrt(U ** 2))

        ax['vmin'] = kwargs.get('vmin', 0)
        ax['vmax'] = kwargs.get('vmax', np.max(max_vector))

        U = U / max_vector

        # Scale factors
        vx_scale = kwargs.get('vx_scale', spacing*self.dx/self.a0)

        # Scaling
        U = vx_scale*U

        if not kwargs['field_is_nan']:
            fig.add_trace(go.Cone(x=X.flatten()/self.a0, 
                                    y=Y.flatten()/self.a0, 
                                    z=Z.flatten()/self.a0, 
                                    u=U.flatten(), 
                                    v=V.flatten(), 
                                    w=W.flatten(), 
                                    colorscale=kwargs['colormap_object'], 
                                    sizemode='scaled', 
                                    sizeref=1, 
                                    showscale=False,
                                    scene=ax['sceneN'],
                                    hovertemplate='<b>x:</b> %{x:.1f} a₀ <br>' +
                                                    '<b>y:</b> %{y:.1f} a₀ <br>' +
                                                    '<b>z:</b> %{z:.1f} a₀ <br>' +
                                                    '<b>ux:</b> %{u:.1f} <br>' +
                                                    '<b>uy:</b> %{v:.1f} <br>' +
                                                    '<b>uz:</b> %{w:.1f} <br>' +
                                                    '<b>|u|:</b> %{customdata:.1f}<extra></extra>',
                                                    customdata=np.sqrt(U**2 + V**2 + W**2).flatten()))

    ###############################################################
    ########### DIMENSION: 3 - VECTOR-DIMENSION: 2 ################
    ###############################################################

    elif self.dim == 3 and vector_field.shape == (2,self.xRes,self.yRes,self.zRes):

        ax = tool_plotly_define_3D_plot_ax(fig, ax)

        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        # Define the vector field
        U = vector_field[0]
        V = vector_field[1]
        W = np.zeros(U.shape)

        # Add spacing
        X,Y,Z,U,V,W = tool_add_spacing_3D(X,Y,Z,U,V,W,spacing)

        # Normalize the vectors
        max_vector = np.max(np.sqrt(U ** 2 + V ** 2))

        ax['vmin'] = kwargs.get('vmin', 0)
        ax['vmax'] = kwargs.get('vmax', np.max(max_vector))

        U = U / max_vector
        V = V / max_vector

        # Scale factors
        vx_scale = kwargs.get('vx_scale', spacing*self.dx/self.a0)
        vy_scale = kwargs.get('vy_scale', spacing*self.dy/self.a0)

        # Scaling
        U = vx_scale*U
        V = vy_scale*V

        if not kwargs['field_is_nan']:
            fig.add_trace(go.Cone(x=X.flatten()/self.a0, 
                                    y=Y.flatten()/self.a0, 
                                    z=Z.flatten()/self.a0, 
                                    u=U.flatten(), 
                                    v=V.flatten(), 
                                    w=W.flatten(), 
                                    colorscale=kwargs['colormap_object'], 
                                    sizemode='scaled', 
                                    sizeref=1, 
                                    showscale=False,
                                    scene=ax['sceneN'],
                                    hovertemplate='<b>x:</b> %{x:.1f} a₀ <br>' +
                                                    '<b>y:</b> %{y:.1f} a₀ <br>' +
                                                    '<b>z:</b> %{z:.1f} a₀ <br>' +
                                                    '<b>ux:</b> %{u:.1f} <br>' +
                                                    '<b>uy:</b> %{v:.1f} <br>' +
                                                    '<b>uz:</b> %{w:.1f} <br>' +
                                                    '<b>|u|:</b> %{customdata:.1f}<extra></extra>',
                                                    customdata=np.sqrt(U**2 + V**2 + W**2).flatten()))

    ###############################################################
    ########### DIMENSION: 3 - VECTOR-DIMENSION: 3 ################
    ###############################################################

    elif self.dim == 3 and vector_field.shape == (3,self.xRes,self.yRes,self.zRes):

        ax = tool_plotly_define_3D_plot_ax(fig, ax)

        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        # Define the vector field
        U = vector_field[0]
        V = vector_field[1]
        W = vector_field[2]

        X,Y,Z,U,V,W = tool_add_spacing_3D(X,Y,Z,U,V,W,spacing)

        max_vector = np.max(np.sqrt(U ** 2 + V ** 2 + W ** 2))
        ax['vmin'] = kwargs.get('vmin', 0)
        ax['vmax'] = kwargs.get('vmax', np.max(max_vector))

        if not kwargs['field_is_nan']:
            fig.add_trace(go.Cone(x=X.flatten()/self.a0, 
                                    y=Y.flatten()/self.a0, 
                                    z=Z.flatten()/self.a0, 
                                    u=U.flatten(), 
                                    v=V.flatten(), 
                                    w=W.flatten(), 
                                    colorscale=kwargs['colormap_object'], 
                                    sizemode='scaled', 
                                    sizeref=1, 
                                    showscale=False,
                                    scene=ax['sceneN'],
                                    hovertemplate='<b>x:</b> %{x:.1f} a₀ <br>' +
                                                    '<b>y:</b> %{y:.1f} a₀ <br>' +
                                                    '<b>z:</b> %{z:.1f} a₀ <br>' +
                                                    '<b>ux:</b> %{u:.1f} <br>' +
                                                    '<b>uy:</b> %{v:.1f} <br>' +
                                                    '<b>uz:</b> %{w:.1f} <br>' +
                                                    '<b>|u|:</b> %{customdata:.1f}<extra></extra>',
                                                    customdata=np.sqrt(U**2 + V**2 + W**2).flatten()))

    ###############################################################
    ###########     NON-VALID DIMENSION            ################
    ###############################################################

    else:
        raise Exception("You have entered an invalid field to the plot_vector_field function.")


    if kwargs['colorbar'] and not(ax['colorbar']) and not(kwargs['field_is_nan']):
        ax['colormap_object'] = kwargs['colormap_object']
        fig.add_trace(tool_plotly_colorbar(ax, type='normal'))
        ax['colorbar'] = True

    kwargs['fig'] = fig
    kwargs['ax'] = ax
    tool_set_plot_axis_properties_plotly(self, **kwargs)
    return fig, ax