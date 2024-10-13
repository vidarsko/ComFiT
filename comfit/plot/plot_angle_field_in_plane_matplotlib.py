from typing import Optional
import numpy as np
import matplotlib

from comfit.plot.plot_complex_field_in_plane_matplotlib import plot_complex_field_in_plane_matplotlib
from comfit.tool.tool_complete_field import tool_complete_field

def plot_angle_field_in_plane_matplotlib(
    self,
    angle_field: np.ndarray,
    normal_vector: Optional[np.ndarray] = None,
    position: Optional[np.ndarray] = None,
    **kwargs
    ) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """Plots the angle field in a plane.

    Args:
        angle_field (numpy.ndarray): The angle field to be plotted.
        normal_vector (array-like, optional): The normal vector of the plane. Default is [0,1,0].
        position (array-like, optional): The position of the plane. Default is the middle of the system.
        **kwargs: Keyword arguments for the plot, see https://comfitlib.com/Plotting/.

    Returns:
        tuple: A tuple containing:
            - matplotlib.figure.Figure: The figure containing the plot.
            - matplotlib.axes.Axes: The axes containing the plot.
    """

    # Check if the vector field is complex
    if np.iscomplexobj(angle_field):
        print("\033[91mWarning: the angle vector field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
        print('Max imaginary part: ', np.max(np.imag(angle_field)))
        angle_field = np.real(angle_field)

    # Extend the field if not a complete array is given
    angle_field = tool_complete_field(self, angle_field)

    complex_field = np.exp(1j * angle_field)

    return plot_complex_field_in_plane_matplotlib(self, complex_field, normal_vector=normal_vector, position=position, **kwargs)
