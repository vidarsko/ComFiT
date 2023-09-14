import numpy as np
from comfit.core.base_system import BaseSystem

class QuantumSystem(BaseSystem):
    def __init__(self, dimension, **kwargs):
        """
        Initializes a quamtum mechanics system evolving according to the Schrödinger equation

        """

        # First initialize the BaseSystem class
        super().__init__(dimension, **kwargs)

        # Type of the system
        self.psi = None
        self.psi_f = None
        self.type = 'BEC'

        # Default simulation parameters
        self.gamma = 0.1  # Dissipation (gamma)
        self.V_ext = 0  # External potential

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)