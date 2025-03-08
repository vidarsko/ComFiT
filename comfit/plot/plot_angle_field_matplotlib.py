import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from comfit.plot.plot_complex_field_matplotlib import plot_complex_field_matplotlib
from comfit.plot import plot_field_matplotlib
from comfit.tool import tool_complete_field



def plot_angle_field_matplotlib(self, 
       angle_field: np.ndarray, 
        **kwargs
        ) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """Plot the angle field.

    A description.

    Args:
        - self (Comfit): a Comfit instance.
        - field (array-like): The angle field values.
        - **kwargs: Additional arguments to pass to the plot_field_matplotlib function.
    
    Returns:
        - tuple: A tuple containing:
            - matplotlib.figure.Figure: The figure containing the plot.
            - matplotlib.axes.Axes: The axes containing the plot.
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