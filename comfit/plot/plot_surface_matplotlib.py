from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np
from skimage.measure import marching_cubes

def plot_surface_matplotlib(
        self: 'BaseSystem',
        **kwargs: Any
        ) -> Axes:
    """Plot the surface of the given field.

    Parameters
    ----------
    self : BaseSystem
        A BaseSystem (or derived) instance.
    \*\*kwargs : Any
        field : ndarray
            3D array containing the field values
        value : float
            Isosurface value
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. Defaults to current axes.
        alpha : float, optional
            Transparency value. Defaults to 0.5.
        color : str, optional
            Surface color. Defaults to 'b'.

    Returns
    -------
    matplotlib.axes.Axes
        The axes containing the surface plot.
    """
    field = kwargs['field']
    value = kwargs['value']
    ax = kwargs.get('ax', plt.gca())
    
    alpha = kwargs.get('alpha', 0.5)
    color = kwargs.get('color', 'b')

    verts, faces, _, _ = marching_cubes(field, value)

    # Extract coordinates
    x = kwargs.get('x', self.x/self.a0).flatten()
    dx = x[1] - x[0]
    xmax = x[-1]+dx

    y = kwargs.get('y', self.y/self.a0).flatten()
    dy = y[1] - y[0]
    ymax = y[-1]+dy

    z = kwargs.get('z', self.z/self.a0).flatten()
    dz = z[1] - z[0]
    zmax = z[-1]+dz

    ax.plot_trisurf(x[0]+verts[:, 0]*dx, 
                    y[0]+verts[:, 1]*dy, 
                    faces, 
                    z[0]+verts[:, 2]*dz, alpha=alpha, color=color)

    return ax