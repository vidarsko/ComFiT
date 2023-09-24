import numpy as np
from comfit.core.base_system import BaseSystem


class PFC(BaseSystem):

    def __init__(self, dimension, **kwargs):
        """
        Nothing here yet
        """

        # First initialize the BaseSystem class
        super().__init__(dimension, **kwargs)

        # Type of the system
        self.type = 'PFC'

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)


class PFCper(PFC):
    def __init__(self,  **kwargs):
        """
        Nothing here yet
        """
        # First initialize the BaseSystem class
        super().__init__(dimension, **kwargs)

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

class PFCtri(PFC):
    def __init__(self, nx, ny, **kwargs):
        """
        Nothing here yet
        """

        # Type of the system
        self.type = 'PFCtri'
        self.dim = 2

        # Default simulation parameters
        self.micro_resolution = [7, 12]
        self.psi0 = -0.3
        self.r = -0.3
        self.dt = 0.1


        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.xRes = nx*self.micro_resolution[0]
        self.yRes = ny*self.micro_resolution[1]

        # First initialize the BaseSystem class
        super().__init__(self.dim, xRes=self.xRes,yRes=self.yRes)

        self.a0 = 2 * np.pi * 2 / np.sqrt(3)
        self.a = [[1 / 2, np.sqrt(3) / 2],[1, 0],[1 / 2, -np.sqrt(3) / 2]]

        # Multiply every entry in self.a with self.a0
        for i in range(len(self.a)):
            for j in range(len(self.a[i])):
                self.a[i][j] *= self.a0

        self.q = [[0, 1], [np.sqrt(3) / 2, -1 / 2], [-np.sqrt(3) / 2, -1 / 2]]





class PFCsq(PFC):
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


class PFCbcc(PFC):
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


class PFCfcc(PFC):
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


class PFCsc(PFC):
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


class PFCper(PFC):
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