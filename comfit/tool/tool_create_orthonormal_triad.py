import numpy as np

def tool_create_orthonormal_triad(t):
    """
    This function creates an orthonormal triad given a vector t.
    The function first normalizes the input vector t, then it creates two orthogonal unit vectors e1 and e2.
    The vectors e1, e2, and the normalized t form an orthonormal triad.

    Args:
    t : numpy array
        The input vector. It does not need to be a unit vector.

    Returns:
    e1, e2, t : tuple of numpy arrays
        The orthonormal triad. Each element of the tuple is a numpy array representing a vector.
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
