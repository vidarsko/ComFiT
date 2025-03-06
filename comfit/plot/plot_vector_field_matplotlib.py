import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from comfit.tool.tool_complete_field import tool_complete_field
from comfit.tool.tool_add_spacing_2D import tool_add_spacing_2D
from comfit.tool.tool_add_spacing_3D import tool_add_spacing_3D

from comfit.tool import tool_matplotlib_define_2D_plot_ax
from comfit.tool import tool_matplotlib_define_3D_plot_ax

from comfit.tool.tool_set_plot_axis_properties_matplotlib import tool_set_plot_axis_properties_matplotlib


def plot_vector_field_matplotlib(self, 
        vector_field, 
        **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
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

    vector_field, fig, ax, kwargs = self.plot_prepare(vector_field, field_type = 'vector', **kwargs)


    ###############################################################
    ########### DIMENSION: 1 - VECTOR-DIMENSION: 1 ################
    ###############################################################
        
    if self.dim == 1 and vector_field.shape == (1,self.xRes):

        ax = tool_matplotlib_define_2D_plot_ax(fig, ax)

        X, Y = np.meshgrid(self.x, np.array([0]), indexing='ij')

        U = np.zeros(X.shape)
        V = np.zeros(X.shape)

        V[:,0] = vector_field[0]

        X,Y,U,V = tool_add_spacing_2D(X,Y,U,V,spacing)

        ax.quiver(X/self.a0, Y/self.a0, U, V, color='blue', angles='xy', scale_units='xy', scale=1)
        kwargs['ylim'] = [np.min(vector_field[0]), np.max(vector_field[0])]
        

    ###############################################################
    ########### DIMENSION: 1 - VECTOR-DIMENSION: 2 ################
    ###############################################################

    elif self.dim == 1 and vector_field.shape == (2,self.xRes):

        ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
        kwargs['plot_is_3D'] = True
            
        X, Y, Z = np.meshgrid(self.x, np.array([0]), np.array([0]), indexing='ij')

        U = np.zeros(X.shape)
        V = np.zeros(X.shape)
        W = np.zeros(X.shape)

        V[:,0,0] = vector_field[0]
        W[:,0,0] = vector_field[1]

        X,Y,Z,U,V,W = tool_add_spacing_3D(X,Y,Z,U,V,W,spacing)

        ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')

        kwargs['ylim'] = [np.min(vector_field[0]), np.max(vector_field[0])]
        delta_y = kwargs['ylim'][1] - kwargs['ylim'][0]

        kwargs['zlim'] = [np.min(vector_field[1]), np.max(vector_field[1])]
        delta_z = kwargs['zlim'][1] - kwargs['zlim'][0]

        if delta_y < 0.15*delta_z:
            kwargs['ylim'] = [kwargs['ylim'][0] - 0.15*delta_z, kwargs['ylim'][1] + 0.15*delta_z]
        
        if delta_z < 0.15*delta_y:
            kwargs['zlim'] = [kwargs['zlim'][0] - 0.15*delta_y, kwargs['zlim'][1] + 0.15*delta_y]

    ###############################################################
    ########### DIMENSION: 1 - VECTOR-DIMENSION: 3 ################
    ###############################################################

    elif self.dim == 1 and vector_field.shape == (3,self.xRes):
        
        ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
        kwargs['plot_is_3D'] = True

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
        vx_scale = kwargs.get('vx_scale', spacing*self.size_x/15/self.a0)
        vy_scale = kwargs.get('vy_scale', 1)
        vz_scale = kwargs.get('vz_scale', 1)

        # Scaling
        U = vx_scale*U
        V = vy_scale*V
        W = vz_scale*W

        ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')

        kwargs['ylim'] = [-1,1]
        kwargs['zlim'] = [-1,1]

    ###############################################################
    ########### DIMENSION: 2 - VECTOR-DIMENSION: 1 ################
    ###############################################################

    elif self.dim == 2 and vector_field.shape == (1,self.xRes,self.yRes):
        
        ax = tool_matplotlib_define_2D_plot_ax(fig, ax)

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

        ax.quiver(X/self.a0, Y/self.a0, U, V, color='blue',
            angles='xy', scale_units='xy', scale=1,
            headwidth=3, headlength=4, headaxislength=3, pivot='middle')

    ###############################################################
    ########### DIMENSION: 2 - VECTOR-DIMENSION: 2 ################
    ###############################################################

    elif self.dim == 2 and vector_field.shape == (2,self.xRes,self.yRes):

        ax = tool_matplotlib_define_2D_plot_ax(fig, ax)

        X, Y = np.meshgrid(self.x, self.y, indexing='ij')

        X, Y, U, V = tool_add_spacing_2D(X,Y,vector_field[0],vector_field[1],spacing)

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


    ###############################################################
    ########### DIMENSION: 2 - VECTOR-DIMENSION: 3 ################
    ###############################################################

    elif self.dim == 2 and vector_field.shape == (3,self.xRes,self.yRes):

        ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
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
        vx_scale = kwargs.get('vx_scale', 2*spacing*self.size_x/max_vector/self.a0)
        vy_scale = kwargs.get('vy_scale', 2*spacing*self.size_y/max_vector/self.a0)
        vz_scale = kwargs.get('vz_scale', spacing/self.a0)

        # Scaling
        U = vx_scale*U
        V = vy_scale*V
        W = vz_scale*W 
        
        ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')
       
        kwargs['axis_equal'] = False
        kwargs['zlim'] = [-spacing,spacing]


    ###############################################################
    ########### DIMENSION: 3 - VECTOR-DIMENSION: 1 ################
    ###############################################################

    elif self.dim == 3 and vector_field.shape == (1,self.xRes,self.yRes,self.zRes):

        ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
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

        ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')
        
    ###############################################################
    ########### DIMENSION: 3 - VECTOR-DIMENSION: 2 ################
    ###############################################################

    elif self.dim == 3 and vector_field.shape == (2,self.xRes,self.yRes,self.zRes):

        ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
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
        

        ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')
        
   
    ###############################################################
    ########### DIMENSION: 3 - VECTOR-DIMENSION: 3 ################
    ###############################################################

    elif self.dim == 3 and vector_field.shape == (3,self.xRes,self.yRes,self.zRes):

        ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
        kwargs['plot_is_3D'] = True

        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        # Define the vector field
        U = vector_field[0]
        V = vector_field[1]
        W = vector_field[2]

        X,Y,Z,U,V,W = tool_add_spacing_3D(X,Y,Z,U,V,W,spacing)



        

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
        

    ###############################################################
    ###########     NON-VALID DIMENSION            ################
    ###############################################################

    else:
        raise Exception("You have entered an invalid field to the plot_vector_field function.")

    kwargs['ax'] = ax
    tool_set_plot_axis_properties_matplotlib(self, **kwargs)
    return fig, ax
