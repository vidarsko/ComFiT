import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib
from typing import Optional
from comfit.tool.tool_set_plot_axis_properties_matplotlib import tool_set_plot_axis_properties_matplotlib


####################################################################################################
########### NOT YET PORTED TO PLOTLY ##############################################################
####################################################################################################

def plot_vector_field_in_plane_plotly(
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

    print("\033[91mWarning: The plot_vector_field_in_plane_plotly function has not been ported to Plotly yet.\033[0m")
    pass

    # # Convert the vector field to a numpy array
    # vector_field = np.array(vector_field)

    # # Extend the field if not a complete array is given
    # # TODO: Make this part more efficient (Vidar 11.03.24)
    # vector_field_copy = []
    # for n in range(vector_field.shape[0]):
    #     vector_field_copy.append(tool_complete_field(self, vector_field[n]))
    # vector_field = np.array(vector_field_copy)

    # # Check if the vector field is complex
    # if np.iscomplexobj(vector_field):
    #     print("\033[91mWarning: the provided vector field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
    #     print('Max imaginary part: ', np.max(np.imag(vector_field)))
    #     vector_field = np.real(vector_field)
    

    # if self.dim != 3:
    #     raise Exception("The plot in plane function is only defined for 3D fields.")

    # # Check if an axis object is provided
    # fig = kwargs.get('fig', plt.gcf())
    # ax = kwargs.get('ax', {'row': 1, 'col': 1, 'nrows': 1, 'ncols': 1})

    # if ax == None:
    #     ax = fig.add_subplot(111, projection='3d')

    # # Default values of position and normal vector
    # if position is None:
    #     position = self.rmid
    
    # if normal_vector is None:
    #     normal_vector = [0,1,0]
    
    # normal_vector = np.array(normal_vector)/np.linalg.norm(normal_vector)
    # height_above_plane = (self.x-position[0])*normal_vector[0] + (self.y-position[1])*normal_vector[1] + (self.z-position[2])*normal_vector[2]

    # verts, faces, _, _ = marching_cubes(height_above_plane, 0)

    # # Assuming field is defined on the same grid as height_above_plane
    # x, y, z = np.mgrid[0:height_above_plane.shape[0], 0:height_above_plane.shape[1], 0:height_above_plane.shape[2]]

    # # Flatten the grid for interpolation
    # points = np.c_[x.ravel(), y.ravel(), z.ravel()]

    # # Extract the vector field
    # if vector_field.shape == (1,self.xRes,self.yRes,self.zRes):
    #     n = 1
    # elif vector_field.shape == (2,self.xRes,self.yRes,self.zRes):
    #     n = 2
    # elif vector_field.shape == (3,self.xRes,self.yRes,self.zRes):
    #     n = 3
    # else:
    #     raise Exception("You have entered an invalid field to the plot_vector_field function.")

    # U = vector_field[0]
    # U_values = U.ravel()
    # U_verts = sp.interpolate.griddata(points, U_values, verts, method='nearest')

    # if n > 1:
    #     V = vector_field[1]
    #     V_values = V.ravel()
    #     V_verts = sp.interpolate.griddata(points, V_values, verts, method='nearest')
    # else:
    #     V_verts = np.zeros(U_verts.shape)
    
    # if n > 2:
    #     W = vector_field[2]
    #     W_values = W.ravel()
    #     W_verts = sp.interpolate.griddata(points, W_values, verts, method='nearest')
    # else:
    #     W_verts = np.zeros(U_verts.shape)

    # # Normalize the vectors
    # max_vector = np.max(np.sqrt(U_verts ** 2 + V_verts ** 2 + W_verts ** 2))
    # U_verts = U_verts / max_vector
    # V_verts = V_verts / max_vector
    # W_verts = W_verts / max_vector

    # # Scale factors
    # vx_scale = kwargs.get('vx_scale', spacing*self.dx)
    # vy_scale = kwargs.get('vy_scale', spacing*self.dy)
    # vz_scale = kwargs.get('vz_scale', spacing*self.dz)

    # # Scaling
    # U_verts = vx_scale*U_verts
    # V_verts = vy_scale*V_verts
    # W_verts = vz_scale*W_verts

    # x = self.xmin+verts[:, 0]*self.dx
    # y = self.ymin+verts[:, 1]*self.dy
    # z = self.zmin+verts[:, 2]*self.dz

    # # Adjust positions based on the coordinate system
    # x = self.xmin + x * self.dx
    # y = self.ymin + y * self.dy
    # z = self.zmin + z * self.dz

    # # Add spacing
    # x = x[::spacing]
    # y = y[::spacing]
    # z = z[::spacing]
    # U_verts = U_verts[::spacing]
    # V_verts = V_verts[::spacing]
    # W_verts = W_verts[::spacing]

    # ax.quiver(x, y, z, U_verts, V_verts, W_verts, color='blue')

    # kwargs['ax'] = ax
    # tool_set_plot_axis_properties_matplotlib(self, **kwargs)
    # return fig, ax