from typing import Tuple

#General packages
import numpy as np

def tool_add_spacing_2D(
        X,
        Y,
        U,
        V,
        spacing
        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Add spacing to a 2D field.

    Parameters
    ----------
    X : np.ndarray
        The x-coordinates.
    Y : np.ndarray
        The y-coordinates.
    U : np.ndarray
        The x-components of the field.
    V : np.ndarray
        The y-components of the field.
    spacing : int
        The spacing to be added.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        The spaced x-coordinates, y-coordinates, x-components of the field, and y-components of the field.
    """

    X = X[::spacing, ::spacing]
    Y = Y[::spacing, ::spacing]
    U = U[::spacing, ::spacing]
    V = V[::spacing, ::spacing]
    return X,Y,U,V