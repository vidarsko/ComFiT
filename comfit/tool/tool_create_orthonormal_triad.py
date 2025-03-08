from typing import Tuple, List, Union

import numpy as np

def tool_create_orthonormal_triad(
        t: Union[np.ndarray, List[float]]
        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create an orthonormal triad given a vector t.
    
    Parameters
    ----------
    t : np.ndarray
        Input vector of shape (3,). Does not need to be normalized.
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        A tuple (e1, e2, t) containing the orthonormal triad vectors.
        Each vector is a normalized numpy array of shape (3,).
    """
    t = np.array(t)

    # Normalize the input vector
    t = t / np.linalg.norm(t)

    # Define two orthogonal vectors
    ey = np.array([0, 1, 0])
    ez = np.array([0, 0, 1])

    # Check the dot product of ey and t
    # If it is less than 0.8, we can use ey to create the orthogonal vector e1
    # Otherwise, we use ez to create e1
    if abs(np.dot(ey, t)) < 0.8:
        e1 = np.cross(ey, t)
    else:
        e1 = np.cross(ez, t)

    # Normalize e1
    e1 /= np.linalg.norm(e1)

    # Create the second orthogonal vector e2 by taking the cross product of t and e1
    e2 = np.cross(t, e1)

    # Return the orthonormal triad
    return e1, e2, t
