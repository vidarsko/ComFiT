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


    def calc_from_amplitudes(self):

        self.psi = self.psi0

        if self.dim == 1:
            X = self.x

        elif self.dim == 2:
            X, Y = np.meshgrid(self.x, self.y, indexing='ij')

        elif self.dim == 3:
            X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        for n in range(self.number_of_reciprocal_modes):

            if self.dim == 1:
                self.psi += 2*self.eta[n]*np.exp(1j*self.q[n][0]*X)

            elif self.dim == 2:
                self.psi += 2*self.eta[n]*np.exp(1j*(self.q[n][0]*X+self.q[n][1]*Y))

            elif self.dim == 3:
                self.psi += 2*self.eta[n]*np.exp(1j*(self.q[n][0]*X+self.q[n][1]*Y+self.q[n][2]*Z))

        self.psi = np.real(self.psi)



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
        self.t = 0
        self.v = 1
        self.dt = 0.1


        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.xRes = nx*self.micro_resolution[0]
        self.yRes = ny*self.micro_resolution[1]


        a0 = 2 * np.pi * 2 / np.sqrt(3)
        self.a = [[1 / 2, np.sqrt(3) / 2],[1, 0],[1 / 2, -np.sqrt(3) / 2]]

        # Multiply every entry in self.a with self.a0
        for i in range(len(self.a)):
            for j in range(len(self.a[i])):
                self.a[i][j] *= a0

        self.q = [[0, 1], [np.sqrt(3) / 2, -1 / 2], [-np.sqrt(3) / 2, -1 / 2]]

        # Set the number of reciprocal modes
        self.number_of_reciprocal_modes = 3
        self.number_of_principal_reciprocal_modes = 3

        # Set the grid
        self.dx = a0/self.micro_resolution[0]
        self.dy = np.sqrt(3)*a0/self.micro_resolution[1]

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes,
                         dx=self.dx, dy=self.dy, dt=self.dt)


        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

        self.A = self.calc_initial_amplitudes()
        self.eta = [self.A, self.A, self.A]

        #Set the elastic constants
        self.el_mu = 3 * self.A ** 2
        self.el_lambda = 3 * self.A ** 2
        self.el_gamma = 0
        self.el_nu = self.el_lambda/((self.dim-1)*self.el_lambda + 2*self.el_mu + self.el_gamma)


    def calc_initial_amplitudes(self):
        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v

        A = (-3 * v * psi0 + np.sqrt(20 * t * v - 15 * v * r + 30 * t * v * psi0 - 36 * v ** 2 * psi0 ** 2)) / (15 * v)
        return A


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