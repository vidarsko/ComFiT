import numpy as np
import scipy as sp

class BaseSystemGet:
    """ Get methods for the base system class"""
    def get_component_from_symmetric_tensor(self, tensor: np.ndarray, i: int, j: int) -> np.ndarray:
        """Gets the i,j component of a symmetric tensor saved in an array structure.
        
        Parameters
        ----------
        tensor : numpy.ndarray
            The symmetric tensor.
        i : int
            The row index.
        j : int
            The column index.
        
        Returns
        -------
        numpy.ndarray
            The i,j component of the tensor.
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

    def get_component_from_symmetric_traceless_tensor(self, tensor: np.ndarray, i: int, j: int) -> np.ndarray:
        """Gets the i,j component of a symmetric traceless tensor saved in an array structure.

        Parameters
        ----------
        tensor : numpy.ndarray
            The symmetric traceless tensor.
        i : int
            The row index.
        j : int
            The column index.
        
        Returns
        -------
        numpy.ndarray
            The i,j component of the tensor.
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

    def get_component_from_antisymmetric_tensor(self, tensor: np.ndarray, i: int, j: int) -> np.ndarray:
        """Gets the i,j component of an antisymmetric tensor saved in an array structure.

        Parameters
        ----------
        tensor : numpy.ndarray
            The antisymmetric tensor, stored as its independent components: a
            length-1 array in 2 dimensions, a length-3 array in 3 dimensions.
        i : int
            The row index.
        j : int
            The column index.

        Returns
        -------
        numpy.ndarray
            The i,j component of the tensor.
        """
        if self.dim == 2:
            if i==j:
                return 0
            return (-1)**i *tensor[0]
        elif self.dim ==3:
            if i ==j:
                return 0
            else:
                return np.sign(j-i)*tensor[i+j-1]
