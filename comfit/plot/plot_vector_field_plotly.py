import numpy as np
import plotly.graph_objects as go
from comfit.tool.tool_set_plot_axis_properties_plotly import tool_set_plot_axis_properties_plotly

from comfit.tool.tool_complete_field import tool_complete_field
from comfit.tool.tool_add_spacing_2D import tool_add_spacing_2D
from comfit.tool.tool_add_spacing_3D import tool_add_spacing_3D

def plot_vector_field_plotly(self, vector_field, **kwargs):
    """Plots a vector field.

    Args:
        vector_field (tuple): Tuple containing the x and y components of the vector field.
        **kwargs: 
            spacing (int, optional): The spacing for the quiver plot. Default is self.xRes//20.
            Keyword arguments for the plot, see https://comfitlib.com/Plotting/.

    Returns:
        Tuple containing
            - The figure containing the plot.
            - The axes containing the plot.
    """
    spacing = kwargs.get('spacing', max(self.xRes//20,1))

    # Extend the field if not a complete array is given
    vector_field_copy = []
    for n in range(len(vector_field)):
        vector_field_copy.append(tool_complete_field(self, vector_field[n]))
    vector_field = np.array(vector_field_copy)
    
    # Keyword arguments
    axis_equal = kwargs.get('axis_equal',True)
    colormap = kwargs.get('colormap', 'viridis')

    # Check if the vector field is complex
    if np.iscomplexobj(vector_field):
        print("\033[91mWarning: the provided vector field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
        print('Max imaginary part: ', np.max(np.imag(vector_field)))
        vector_field = np.real(vector_field)
     
    kwargs['plot_is_3D'] = False

    ###############################################################
    ########### DIMENSION: 1 - VECTOR-DIMENSION: 1 ################
    ###############################################################
        
    if self.dim == 1 and vector_field.shape == (1,self.xRes):

        X, Y = np.meshgrid(self.x, np.array([0]), indexing='ij')

        U = np.zeros(X.shape)
        V = np.zeros(X.shape)

        V[:,0] = vector_field[0]

        X,Y,U,V = tool_add_spacing_2D(X,Y,U,V,spacing)

            
        fig = kwargs.get('fig', go.Figure())
        fig = self.plot_field(vector_field[0], **kwargs)

        

    ###############################################################
    ########### DIMENSION: 1 - VECTOR-DIMENSION: 2 ################
    ###############################################################

    elif self.dim == 1 and vector_field.shape == (2,self.xRes):

        kwargs['plot_is_3D'] = True
            
        X, Y, Z = np.meshgrid(self.x, np.array([0]), np.array([0]), indexing='ij')

        U = np.zeros(X.shape)
        V = np.zeros(X.shape)
        W = np.zeros(X.shape)

        V[:,0,0] = vector_field[0]
        W[:,0,0] = vector_field[1]

        X,Y,Z,U,V,W = tool_add_spacing_3D(X,Y,Z,U,V,W,spacing)

        fig = kwargs.get('fig', go.Figure())
        fig.add_trace(go.Cone(x=X.flatten()/self.a0, 
                                y=Y.flatten()/self.a0, 
                                z=Z.flatten()/self.a0, 
                                u=U.flatten(), 
                                v=V.flatten(), 
                                w=W.flatten(), 
                                colorscale='Viridis', 
                                sizemode='scaled', 
                                sizeref=1, 
                                showscale=True,
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
        
        kwargs['plot_is_3D'] = True


        fig = kwargs.get('fig', go.Figure())
        
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

        fig = kwargs.get('fig', go.Figure())
        fig.add_trace(go.Cone(x=X.flatten()/self.a0, 
                        y=Y.flatten()/self.a0, 
                        z=Z.flatten()/self.a0, 
                        u=U.flatten(), 
                        v=V.flatten(), 
                        w=W.flatten(), 
                        colorscale='Viridis', 
                        sizemode='scaled', 
                        sizeref=1, 
                        showscale=True, 
                        hovertemplate='<b>x:</b> %{x:.1f} a₀ <br>' +
                                      '<b>ux:</b> %{u:.1f} <br>' +
                                      '<b>uy:</b> %{v:.1f} <br>' +
                                      '<b>uz:</b> %{w:.1f} <br>' +
                                      '<b>|u|:</b> %{customdata:.1f}<extra></extra>',
                                      customdata=np.sqrt(U**2 + V**2 + W**2).flatten()))

        kwargs['axis_equal'] = False
        # kwargs['ylim'] = [-1,1]
        # kwargs['zlim'] = [-1,1]

    ###############################################################
    ########### DIMENSION: 2 - VECTOR-DIMENSION: 1 ################
    ###############################################################

    elif self.dim == 2 and vector_field.shape == (1,self.xRes,self.yRes):
  

        X, Y = np.meshgrid(self.x, self.y, indexing='ij')

        U = vector_field[0]
        V = np.zeros(X.shape)

        X,Y,U,V = tool_add_spacing_2D(X,Y,U,V,spacing)

        # Normalize the vectors
        max_vector = np.max(np.sqrt(U ** 2))
        
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

        colorbar = kwargs.get('colorbar', True)
        
        fig.add_trace(go.Scatter(
        x=X.flatten()/self.a0,
        y=Y.flatten()/self.a0,
        mode='markers',
        marker=dict(symbol='arrow', 
            angle=90-angle.flatten()*180/np.pi, 
            size=2*spacing*magnitude_normalized.flatten(), 
            sizemode='diameter',
            color=magnitude.flatten(), 
            colorscale='Viridis', 
            showscale=colorbar,
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

        X, Y = np.meshgrid(self.x, self.y, indexing='ij')

        X, Y, U, V = tool_add_spacing_2D(X,Y,vector_field[0],vector_field[1],spacing)
        
        fig = kwargs.get('fig', go.Figure())

        u = U.flatten()
        v = V.flatten()
        
        magnitude = np.sqrt(u**2 + v**2)
        magnitude_max = max(np.max(magnitude),1e-12)
        magnitude_normalized = magnitude/magnitude_max

        angle = np.arctan2(v, u)
        direction = np.array([np.cos(angle), np.sin(angle)]).T

        colorbar = kwargs.get('colorbar', True)
        
        fig.add_trace(go.Scatter(
        x=X.flatten()/self.a0,
        y=Y.flatten()/self.a0,
        mode='markers',
        marker=dict(symbol='arrow', 
            angle=90-angle.flatten()*180/np.pi, 
            size=2*spacing*magnitude_normalized.flatten(), 
            sizemode='diameter',
            color=magnitude.flatten(), 
            colorscale='Viridis', 
            showscale=colorbar,
            cmin=0,
            cmax=magnitude_max,
            line=dict(color='black')
            ),
            hovertemplate='<b>x:</b> %{x:.2f}a0<br>' +
                        '<b>y:</b> %{y:.2f}a0<br>' +
                        '<b>ux:</b> %{customdata[0]:.2f}<br>' +
                        '<b>uy:</b> %{customdata[1]:.2f}<br>' +
                        '<b>|u|:</b> %{customdata[2]:.2f}<extra></extra>',
            customdata=np.stack((u.flatten(), v.flatten(), magnitude.flatten()), axis=-1)  # Adding ux, uy and u as customdata
        )
        )


    ###############################################################
    ########### DIMENSION: 2 - VECTOR-DIMENSION: 3 ################
    ###############################################################

    elif self.dim == 2 and vector_field.shape == (3,self.xRes,self.yRes):

        kwargs['plot_is_3D'] = True
        
        X, Y, Z = np.meshgrid(self.x, self.y, np.array([0]), indexing='ij')
        U = np.zeros(X.shape)
        V = np.zeros(X.shape)
        W = np.zeros(X.shape)
        
        U[:,:,0] = vector_field[0]
        V[:,:,0] = vector_field[1]
        W[:,:,0] = vector_field[2]

        X,Y,Z,U,V,W = tool_add_spacing_3D(X,Y,Z,U,V,W,spacing)

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

        fig = kwargs.get('fig', go.Figure())
        fig.add_trace(go.Cone(x=X.flatten()/self.a0, 
                            y=Y.flatten()/self.a0, 
                            z=Z.flatten()/self.a0, 
                            u=U.flatten(), 
                            v=V.flatten(), 
                            w=W.flatten(), 
                            colorscale='Viridis', 
                            sizemode='scaled', 
                            sizeref=1, 
                            showscale=True,
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

        kwargs['plot_is_3D'] = True

        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')              

        # Define the vector field
        U = vector_field[0]
        V = np.zeros(U.shape)
        W = np.zeros(U.shape)

        # Add spacing
        X,Y,Z,U,V,W = tool_add_spacing_3D(X,Y,Z,U,V,W,spacing)

        # Normalize the vectors
        max_vector = np.max(np.sqrt(U ** 2))
        U = U / max_vector

        # Scale factors
        vx_scale = kwargs.get('vx_scale', spacing*self.dx/self.a0)

        # Scaling
        U = vx_scale*U

        fig = kwargs.get('fig', go.Figure())
        fig.add_trace(go.Cone(x=X.flatten()/self.a0, 
                                y=Y.flatten()/self.a0, 
                                z=Z.flatten()/self.a0, 
                                u=U.flatten(), 
                                v=V.flatten(), 
                                w=W.flatten(), 
                                colorscale='Viridis', 
                                sizemode='scaled', 
                                sizeref=1, 
                                showscale=True,
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

        kwargs['plot_is_3D'] = True

        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        # Define the vector field
        U = vector_field[0]
        V = vector_field[1]
        W = np.zeros(U.shape)

        # Add spacing
        X,Y,Z,U,V,W = tool_add_spacing_3D(X,Y,Z,U,V,W,spacing)

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


        fig = kwargs.get('fig', go.Figure())
        fig.add_trace(go.Cone(x=X.flatten()/self.a0, 
                                y=Y.flatten()/self.a0, 
                                z=Z.flatten()/self.a0, 
                                u=U.flatten(), 
                                v=V.flatten(), 
                                w=W.flatten(), 
                                colorscale='Viridis', 
                                sizemode='scaled', 
                                sizeref=1, 
                                showscale=True,
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

        kwargs['plot_is_3D'] = True

        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        # Define the vector field
        U = vector_field[0]
        V = vector_field[1]
        W = vector_field[2]

        X,Y,Z,U,V,W = tool_add_spacing_3D(X,Y,Z,U,V,W,spacing)

        fig = kwargs.get('fig', go.Figure())
        fig.add_trace(go.Cone(x=X.flatten()/self.a0, 
                                y=Y.flatten()/self.a0, 
                                z=Z.flatten()/self.a0, 
                                u=U.flatten(), 
                                v=V.flatten(), 
                                w=W.flatten(), 
                                colorscale='Viridis', 
                                sizemode='scaled', 
                                sizeref=1, 
                                showscale=True,
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



    kwargs['fig'] = fig
    tool_set_plot_axis_properties_plotly(self, **kwargs)
    return fig