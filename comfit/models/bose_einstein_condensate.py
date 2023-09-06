import numpy as np
from comfit.core.base_system import BaseSystem

class BEC(BaseSystem):
    def __init__(self, dimension, x_resolution, **kwargs):
        """
        Initializes a system to simulate a Bose-Einstein Condensate using the Gross-Pitaevskii equation.

        Parameters:
        - dimension : int
            The dimension of the system.
        - x_resolution : int
            The resolution along the x-axis.
        - kwargs : dict, optional
            Optional keyword arguments to set additional parameters.

        Returns:
        - BEC object
            The system object representing the BEC simulation.

        Example:
        bec = BEC(3, 100, dissipation_factor=0.5)
        Creates a BEC system with 3 dimensions and an x-resolution of 100. The dissipation factor is set to 0.5.
        """

        # First initialize the BaseSystem class
        super().__init__(dimension, x_resolution)

        # Type of the system
        self.type = 'BEC'

        # Default simulation parameters
        self.gamma = 0        # Dissipation (gamma)
        self.V_ext = 0        # External potential

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)
