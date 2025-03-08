# Typing imports
from typing import TYPE_CHECKING, Any, Tuple
if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

# General packages
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Comfit packages
from comfit.plot import plot_field_matplotlib, plot_complex_field_matplotlib
from comfit.tool import tool_complete_field

def plot_angle_field_matplotlib(
        angle_field: np.ndarray,
        **kwargs: Any
        ) -> Tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """Plot an angle field using matplotlib.

    Parameters
    ----------
    angle_field : np.ndarray
        The angle field values in radians.
    \*\*kwargs : Any
        Keyword arguments for the plot, see https://comfitlib.com/Plotting/.

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        The figure and axes objects containing the plot.
    """

    # Normalize around 0
    angle_field = np.mod(angle_field + np.pi, 2 * np.pi) - np.pi        

    if self.dim == 1:
        if 'vlim' in kwargs:
            vlim = kwargs['vlim']
        else:
            kwargs['vlim'] = [-np.pi, np.pi]
            kwargs['yticks'] = [-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi]
            kwargs['yticklabels'] = [r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$']

        
        return plot_field_matplotlib(self, angle_field, **kwargs)
    
    elif self.dim > 1:
        complex_field = np.exp(1j * angle_field)

        kwargs['plot_method'] = 'phase_angle'

        return plot_complex_field_matplotlib(self, complex_field, **kwargs)