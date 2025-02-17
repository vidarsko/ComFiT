import numpy as np
from comfit.tool import tool_complete_field
from comfit.plot import plot_field_plotly, plot_complex_field_plotly

def plot_angle_field_plotly(self, angle_field: np.ndarray, **kwargs):
    """Plot the angle field.

    Args:
        field (array-like): The angle field values.
        ax (matplotlib.axes.Axes, optional): The axes to plot the angle field on. If not provided, a new subplot will be created.
    
    Returns:
        The axes containing the plot. (matplotlib.axes.Axes)
    """

    # # Check if the vector field is complex
    if np.iscomplexobj(angle_field):
        print("\033[91mWarning: the angle vector field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
        print('Max imaginary part: ', np.max(np.imag(angle_field)))
        angle_field = np.real(angle_field)

    # Extend the field if not a complete array is given
    angle_field = tool_complete_field(self, angle_field)

    # Normalize around 0
    angle_field = np.mod(angle_field + np.pi, 2 * np.pi) - np.pi        

    if self.dim == 1:
        if 'vlim' in kwargs:
            vlim = kwargs['vlim']
        else:
            kwargs['vlim'] = [-np.pi, np.pi]
            kwargs['yticks'] = [-np.pi, -2 * np.pi / 3, -np.pi / 3, 0, np.pi / 3, 2 * np.pi / 3, np.pi]
            kwargs['yticklabels'] = [r'$-\pi$', r'$-2\pi/3$', r'$-\pi/3$', r'$0$', r'$\pi/3$', r'$2\pi/3$', r'$\pi$']

        
        return plot_field_plotly(self, angle_field, **kwargs)
    
    elif self.dim > 1:
        complex_field = np.exp(1j * angle_field)

        kwargs['plot_method'] = 'phase_angle'

        return plot_complex_field_plotly(self, complex_field, **kwargs)