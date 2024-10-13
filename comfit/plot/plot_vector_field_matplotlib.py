
from comfit.tool.tool_complete_field import tool_complete_field
from comfit.tool.tool_add_spacing_2D import tool_add_spacing_2D
from comfit.tool.tool_add_spacing_3D import tool_add_spacing_3D


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


    # Check if an axis object is provided
    fig = kwargs.get('fig', plt.gcf())
    ax = kwargs.get('ax', None)          


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


        if ax == None:
            fig.clf()
            ax = plt.gcf().add_subplot(111)

        ax.quiver(X/self.a0, Y/self.a0, U, V, color='blue', angles='xy', scale_units='xy', scale=1)
        kwargs['ylim'] = [np.min(vector_field[0]), np.max(vector_field[0])]
        

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


        if ax == None:
            fig.clf()
            ax = fig.add_subplot(111, projection='3d')

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
        
        kwargs['plot_is_3D'] = True

       
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

        ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')

        kwargs['ylim'] = [-1,1]
        kwargs['zlim'] = [-1,1]

    ###############################################################
    ########### DIMENSION: 2 - VECTOR-DIMENSION: 1 ################
    ###############################################################

    elif self.dim == 2 and vector_field.shape == (1,self.xRes,self.yRes):
        
      
        if ax == None:
            fig.clf()
            ax = plt.gcf().add_subplot(111)

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

        if ax == None:
            fig.clf()
            ax = plt.gcf().add_subplot(111)

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

        
        if ax == None:
            fig.clf()
            ax = fig.add_subplot(111, projection='3d')
        ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')
       
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


        if ax == None:
            fig.clf()
            ax = fig.add_subplot(111, projection='3d')
        ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')
        
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

      
        if ax == None:
            fig.clf()
            ax = fig.add_subplot(111, projection='3d')

        ax.quiver(X/self.a0, Y/self.a0, Z/self.a0, U, V, W, color='blue')
        
   
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



        if ax == None:
            fig.clf()
            ax = fig.add_subplot(111, projection='3d')

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
    plot_tool_set_axis_properties_matplotlib(self, **kwargs)
    return fig, ax
