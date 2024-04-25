import numpy as np
import scipy as sp

class BaseSystemGet:
    # GET FUNCTIONS

    def get_sym(self,tensor,i,j):
        """
        Gets the i,j component of a symmetric tensor saved in an array structure.

        Args: 
            tensor (numpy.ndarray): The symmetric tensor.
            i (int): The row index.
            j (int): The column index.
        
        Output:
            (numpy.ndarray) The i,j component of the tensor.
        """

        if self.dim == 2:
            if i == 0:
                return tensor[0] if j == 0 else tensor[1]
            elif i == 1:
                return tensor[1] if j == 0 else tensor[2]

        elif self.dim == 3:
            if i == 0:
                return tensor[0] if j == 0 else tensor[1] if j == 1 else tensor[2]
            elif i == 1:
                return tensor[1] if j == 0 else tensor[3] if j == 1 else tensor[4]
            elif i == 2:
                return tensor[2] if j == 0 else tensor[4] if j == 1 else tensor[5]

    def get_sym_tl(self,tensor,i,j):
        """
        Gets the i,j component of a symmetric traceless tensor saved in an array structure.

        Args:
            tensor (numpy.ndarray): The symmetric traceless tensor.
            i (int): The row index.
            j (int): The column index.
        
        Output:
            (numpy.ndarray) The i,j component of the tensor.
        """
        if self.dim == 2:
            if i == 0:
                return tensor[0] if j == 0 else tensor[1]
            elif i == 1:
                return tensor[1] if j == 0 else -tensor[0]
            
        elif self.dim == 3:
            if i == 0:
                return tensor[0] if j == 0 else tensor[1] if j == 1 else tensor[2]
            elif i == 1:
                return tensor[1] if j == 0 else tensor[3] if j == 1 else tensor[4]
            elif i == 2:
                return tensor[2] if j == 0 else tensor[4] if j == 1 else -(tensor[0] + tensor[3])

    def get_anti_sym(self,omega,i,j):
        """
        Gets the i,j component of an anti-symmetric tensor saved in an array structure.

        Args:
            omega (numpy.ndarray): The anti-symmetric tensor.
            i (int): The row index.
            j (int): The column index.
        
        Output:
            (numpy.ndarray) The i,j component of the tensor.
        """
        # TODO: I don't like that the input vector is a scalar field in 2 dimensions. (Vidar 11.03.24)
        if self.dim == 2:
            if i==j:
                return 0
            return (-1)**i *omega
        elif self.dim ==3:
            if i ==j:
                return 0
            else:
                return np.sign(j-i)*omega[i+j-1]

