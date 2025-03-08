from typing import Tuple

#General packages
import numpy as np

def tool_add_spacing_3D(
        X : np.ndarray,
        Y : np.ndarray,
        Z : np.ndarray,
        U : np.ndarray,
        V : np.ndarray,
        W : np.ndarray,
        spacing : int
        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Add spacing to a 3D field.

    When working with plot vectors, it is often useful to reduce the number of 
    vectors to make the plot more readable. This function allows you to do that 
    by adding spacing to the field.

    Parameters
    ----------
    X : np.ndarray
        The x-coordinates.
    Y : np.ndarray
        The y-coordinates.
    Z : np.ndarray
        The z-coordinates.
    U : np.ndarray
        The x-components of the field.
    V : np.ndarray
        The y-components of the field.
    W : np.ndarray
        The z-components of the field.
    spacing : int
        The spacing to be added.
    
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        The spaced x-, y-, z-coordinates, x-, y-, and z-components of the field.
    """
    X = X[::spacing, ::spacing, ::spacing]
    Y = Y[::spacing, ::spacing, ::spacing]
    Z = Z[::spacing, ::spacing, ::spacing]
    U = U[::spacing, ::spacing, ::spacing]
    V = V[::spacing, ::spacing, ::spacing]
    W = W[::spacing, ::spacing, ::spacing]
    return X,Y,Z,U,V,W  