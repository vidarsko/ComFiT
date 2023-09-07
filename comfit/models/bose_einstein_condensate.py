import numpy as np
from comfit.core.base_system import BaseSystem


class BEC(BaseSystem):
    def __init__(self, dimension, **kwargs):
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
        super().__init__(dimension, **kwargs)

        # Type of the system
        self.type = 'BEC'

        # Default simulation parameters
        self.gamma = 0  # Dissipation (gamma)
        self.V_ext = 0  # External potential

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

    # Initial conditions
    def set_initial_condition_disordered(self, noise_strength=0.01):
        """
        Sets disordered initial condition for the BEC with some thermal flcutiations

        :param noise_strength:
        :return:
        """

        if self.dim==1:
            self.psi = np.random.rand(self.xRes) - 0.5 + \
                       1j * np.random.rand(self.xRes) - 0.5j

        elif self.dim==2:
            self.psi = np.random.rand(self.xRes, self.yRes)-0.5 + \
                       1j*np.random.rand(self.xRes, self.yRes)-0.5j

        elif self.dim==3:
            self.psi = np.random.rand(self.xRes, self.yRes,self.zRes) - 0.5 + \
                       1j * np.random.rand(self.xRes, self.yRes,self.zRes) - 0.5j
        else:
            raise Exception("Code for this dimension has not yet been implemented.")

        self.psi = noise_strength*self.psi