import numpy as np
from comfit.core.base_system import BaseSystem


class PFC(BaseSystem):

    def __init__(self, dimension, x_resolution, **kwargs):
        """
        Nothing here yet
        """

        # First initialize the BaseSystem class
        super().__init__(dimension, x_resolution)

        # Type of the system
        self.type = 'PFC'

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)
