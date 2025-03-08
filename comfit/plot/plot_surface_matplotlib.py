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
    x = (self.xmin+verts[:, 0]*self.dx)/self.a0
    y = (self.ymin+verts[:, 1]*self.dy)/self.a0
    z = (self.zmin+verts[:, 2]*self.dz)/self.a0
    ax.plot_trisurf(x, y, faces, z, alpha=alpha, color=color)

    return ax