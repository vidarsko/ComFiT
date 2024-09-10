import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from skimage.measure import marching_cubes

def plot_surface_matplotlib(self, **kwargs) -> matplotlib.axes.Axes:
    """Plots the surface of the given field.

    Args:
        **kwargs: Keyword arguments for the plot.
    
    Returns:
        The axes containing the plot.
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