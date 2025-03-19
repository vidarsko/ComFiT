from typing import TYPE_CHECKING, Any, Tuple, Union, List

if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

# Standard library imports
import numpy as np

# Third-party library imports
import matplotlib.pyplot as plt
import matplotlib

# Local application imports
from comfit.tool import (
    tool_complete_field,
    tool_add_spacing_2D,
    tool_add_spacing_3D,
    tool_matplotlib_define_2D_plot_ax,
    tool_matplotlib_define_3D_plot_ax, 
    tool_set_plot_axis_properties_matplotlib
)

def plot_vector_field_matplotlib(
        self: 'BaseSystem',
        vector_field: Union[np.ndarray, List],
        **kwargs: Any
        ) -> Tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """Plot the vector field using matplotlib.

    Parameters
    ----------
    self : BaseSystem
        A BaseSystem (or derived) instance.
    vector_field : np.ndarray
        Array containing the components of the vector field
    kwargs : Any
        spacing : int, optional. The spacing for the quiver plot which defaults 
        to self.xRes//20. Additional keyword arguments for plot customization, 
        see https://comfitlib.com/Plotting/

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        Figure and axes containing the plot

    Raises
    ------
    Exception
        If an invalid field is provided to the plot_vector_field function
    """

    spacing = kwargs.get('spacing', max(self.xRes//20,1))

    vector_field, fig, ax, kwargs = self.plot_prepare(vector_field, field_type = 'vector', **kwargs)

    # Extract coordinates
    x = kwargs.get('x', self.x/self.a0).flatten()
    dx = x[1] - x[0]
    xmin = x[0]
    xmax = x[-1]+dx
    
    if self.dim > 1:
        y = kwargs.get('y', self.y/self.a0).flatten()
        dy = y[1] - y[0]
        ymin = y[0]
        ymax = y[-1]+dy

    if self.dim > 2:
        z = kwargs.get('z', self.z/self.a0).flatten()
        dz = z[1] - z[0]
        zmin = z[0]
        zmax = z[-1]+dz


    ###############################################################
    ########### DIMENSION: 1 - VECTOR-DIMENSION: 1 ################
    ###############################################################
        
    if self.dim == 1 and vector_field.shape == (1,self.xRes):

        ax = tool_matplotlib_define_2D_plot_ax(fig, ax)

        X, Y = np.meshgrid(x, np.array([0]), indexing='ij')

        U = np.zeros(X.shape)
        V = np.zeros(X.shape)

        V[:,0] = vector_field[0]

        X,Y,U,V = tool_add_spacing_2D(X,Y,U,V,spacing)

        ax.quiver(X, Y, U, V, color='blue', angles='xy', scale_units='xy', scale=1)
        kwargs['ylim'] = [np.min(vector_field[0]), np.max(vector_field[0])]
        

    ###############################################################
    ########### DIMENSION: 1 - VECTOR-DIMENSION: 2 ################
    ###############################################################

    elif self.dim == 1 and vector_field.shape == (2,self.xRes):

        ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
        kwargs['plot_is_3D'] = True
            
        X, Y, Z = np.meshgrid(x, np.array([0]), np.array([0]), indexing='ij')

        U = np.zeros(X.shape)
        V = np.zeros(X.shape)
        W = np.zeros(X.shape)

        V[:,0,0] = vector_field[0]
        W[:,0,0] = vector_field[1]

        X,Y,Z,U,V,W = tool_add_spacing_3D(X,Y,Z,U,V,W,spacing)

        ax.quiver(X, Y, Z, U, V, W, color='blue')

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

        X, Y, Z = np.meshgrid(x, np.array([0]), np.array([0]), indexing='ij')

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

        ax.quiver(X, Y, Z, U, V, W, color='blue')

        kwargs['ylim'] = [-1,1]
        kwargs['zlim'] = [-1,1]

    ###############################################################
    ########### DIMENSION: 2 - VECTOR-DIMENSION: 1 ################
    ###############################################################

    elif self.dim == 2 and vector_field.shape == (1,self.xRes,self.yRes):
        
        ax = tool_matplotlib_define_2D_plot_ax(fig, ax)

        X, Y = np.meshgrid(x, y, indexing='ij')

        U = vector_field[0]
        V = np.zeros(X.shape)

        X,Y,U,V = tool_add_spacing_2D(X,Y,U,V,spacing)

        # Normalize the vectors
        max_vector = np.max(np.sqrt(U ** 2))
        
        U = U / max_vector

        # Scale factors
        vx_scale = kwargs.get('vx_scale', spacing*dx)

        # Scaling
        U = vx_scale*U

        ax.quiver(X, Y, U, V, color='blue',
            angles='xy', scale_units='xy', scale=1,
            headwidth=3, headlength=4, headaxislength=3, pivot='middle')

    ###############################################################
    ########### DIMENSION: 2 - VECTOR-DIMENSION: 2 ################
    ###############################################################

    elif self.dim == 2 and vector_field.shape == (2,self.xRes,self.yRes):

        ax = tool_matplotlib_define_2D_plot_ax(fig, ax)

        X, Y = np.meshgrid(x, y, indexing='ij')

        X, Y, U, V = tool_add_spacing_2D(X,Y,vector_field[0],vector_field[1],spacing)

        max_vector = np.max(np.sqrt(U ** 2 + V ** 2))
    
        U = U / max_vector
        V = V / max_vector

        # Scale factors
        vx_scale = kwargs.get('vx_scale', 1.0*spacing*dx)
        vy_scale = kwargs.get('vy_scale', 1.0*spacing*dy)

        # Scaling
        U = vx_scale*U
        V = vy_scale*V

        ax.quiver(X, Y, U, V, color='blue', 
            angles='xy', scale_units='xy', scale=1,
            headwidth=3, headlength=4, headaxislength=3, pivot='middle')


    ###############################################################
    ########### DIMENSION: 2 - VECTOR-DIMENSION: 3 ################
    ###############################################################

    elif self.dim == 2 and vector_field.shape == (3,self.xRes,self.yRes):

        ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
        kwargs['plot_is_3D'] = True
        
        X, Y, Z = np.meshgrid(x, y, np.array([0]), indexing='ij')
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
        vx_scale = kwargs.get('vx_scale', 2*spacing*self.size_x/max_vector)
        vy_scale = kwargs.get('vy_scale', 2*spacing*self.size_y/max_vector)
        vz_scale = kwargs.get('vz_scale', spacing)

        # Scaling
        U = vx_scale*U
        V = vy_scale*V
        W = vz_scale*W 
        
        ax.quiver(X, Y, Z, U, V, W, color='blue')
       
        kwargs['axis_equal'] = False
        kwargs['zlim'] = [-spacing,spacing]


    ###############################################################
    ########### DIMENSION: 3 - VECTOR-DIMENSION: 1 ################
    ###############################################################

    elif self.dim == 3 and vector_field.shape == (1,self.xRes,self.yRes,self.zRes):

        ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
        kwargs['plot_is_3D'] = True

        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')              

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
        vx_scale = kwargs.get('vx_scale', spacing*dx)

        # Scaling
        U = vx_scale*U

        ax.quiver(X, Y, Z, U, V, W, color='blue')
        
    ###############################################################
    ########### DIMENSION: 3 - VECTOR-DIMENSION: 2 ################
    ###############################################################

    elif self.dim == 3 and vector_field.shape == (2,self.xRes,self.yRes,self.zRes):

        ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
        kwargs['plot_is_3D'] = True

        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

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
        vx_scale = kwargs.get('vx_scale', spacing*dx)
        vy_scale = kwargs.get('vy_scale', spacing*dy)

        # Scaling
        U = vx_scale*U
        V = vy_scale*V
        

        ax.quiver(X, Y, Z, U, V, W, color='blue')
        
   
    ###############################################################
    ########### DIMENSION: 3 - VECTOR-DIMENSION: 3 ################
    ###############################################################

    elif self.dim == 3 and vector_field.shape == (3,self.xRes,self.yRes,self.zRes):

        ax = tool_matplotlib_define_3D_plot_ax(fig, ax)
        kwargs['plot_is_3D'] = True

        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

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
        vx_scale = kwargs.get('vx_scale', spacing*dx)
        vy_scale = kwargs.get('vy_scale', spacing*dy)
        vz_scale = kwargs.get('vz_scale', spacing*dz)

        # Scaling
        U = vx_scale*U
        V = vy_scale*V
        W = vz_scale*W
        
        ax.quiver(X, Y, Z, U, V, W, color='blue')
        

    ###############################################################
    ###########     NON-VALID DIMENSION            ################
    ###############################################################

    else:
        raise Exception("You have entered an invalid field to the plot_vector_field function.")

    kwargs['ax'] = ax
    tool_set_plot_axis_properties_matplotlib(self, **kwargs)
    return fig, ax
