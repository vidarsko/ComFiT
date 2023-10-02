import numpy as np
from comfit.core.base_system import BaseSystem


class PFC(BaseSystem):

    def __init__(self, dimension, **kwargs):
        """
        Nothing here yet
        """

        # First initialize the BaseSystem class
        super().__init__(dimension, **kwargs)

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.dislocation_charges = np.array(
            [[np.round(np.dot(an, qn) / (2 * np.pi), decimals=8) for qn in self.q] for an in self.a])


    #Calculation methods
    def calc_from_amplitudes(self,eta=None):

        self.psi = self.psi0

        if eta is None:
            eta = self.eta0

        if self.dim == 1:
            X = self.x

        elif self.dim == 2:
            X, Y = np.meshgrid(self.x, self.y, indexing='ij')

        elif self.dim == 3:
            X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        for n in range(self.number_of_reciprocal_modes):

            if self.dim == 1:
                self.psi += 2*eta[n]*np.exp(1j*self.q[n][0]*X)

            elif self.dim == 2:
                self.psi += 2*eta[n]*np.exp(1j*(self.q[n][0]*X+self.q[n][1]*Y))

            elif self.dim == 3:
                self.psi += 2*eta[n]*np.exp(1j*(self.q[n][0]*X+self.q[n][1]*Y+self.q[n][2]*Z))

        self.psi = np.real(self.psi)
        self.psi_f = np.fft.fftn(self.psi)



    def calc_evolution_integrating_factors_PFC_f(self):
        k2 = self.calc_k2()

        if self.type == 'PFCtri':
            omega_f = self.r + (1 - k2)**2
        else:
            raise Exception("Not yet configured.")

        integrating_factors_f = [0,0,0]

        integrating_factors_f[0] = np.exp(-k2 * omega_f * self.dt)
        If1 = integrating_factors_f[0]
        integrating_factors_f[1] = (If1 - 1) / omega_f
        with np.errstate(divide='ignore',invalid='ignore'):
            #This calculation is not valid for the k=0 point, but this will be set explicitly afterwards.
            #That is why it is ignored
            integrating_factors_f[2] = 1. / (-k2 * self.dt * omega_f**2) * (If1 - 1 + k2 * omega_f * self.dt)



        return integrating_factors_f

    def calc_non_linear_evolution_term_PFC_f(self, psi):
        return np.fft.fftn(self.t*psi**2 + self.v*psi**3)


    # Evolution functions
    def evolve_PFC(self,number_of_steps):
        integrating_factors_f = self.calc_evolution_integrating_factors_PFC_f()

        for n in range(number_of_steps):
            self.psi, self.psi_f = self.evolve_ETDRK2_loop(integrating_factors_f,
                                                           self.calc_non_linear_evolution_term_PFC_f,
                                                           self.psi, self.psi_f)
            self.psi = np.real(self.psi)
            self.psi_f = np.fft.fftn(self.psi)


    #Initial configuration methods
    def conf_insert_dislocation(self,eta=None,x=None,y=None,dislocation_type=1):

        if not(self.dim==2):
            raise Exception("The dimension of the system must be 2 for a single point dislocation.")

        if x==None:
            x = self.xmid
        if y==None:
            y = self.ymid
        if eta==None:
            eta = self.eta0

        sn = self.dislocation_charges[dislocation_type-1]
        for n in range(self.number_of_reciprocal_modes):
            if sn[n] != 0:
                eta[n] *= self.calc_angle_field_single_vortex([x,y],charge=sn[n])

        return eta

    def conf_insert_dislocation_dipole(self,eta=None,x1=None,y1=None,x2=None,y2=None,dislocation_type=1):
        """
        Inserts a  dislocation dipole in the system corresponding to dislocation type and its negative

        :param eta: The amplitudes currently being manipulated (default: self.eta0)
        :param x1: The x-position of the first dislocation (default: 1/3 xmax)
        :param y1: The y-position of the first dislocation (default: 1/2 ymax)
        :param x2: The x-position of the second dislocation (default: 2/3 xmax)
        :param y2: The y-position of the second dislocation (default: 1/2 ymax)
        :param dislocation_type: the type of dislocation dipole (default: 1)
        :return:
        """
        #TODO: Improve code documentation

        if not(self.dim==2):
            raise Exception("The dimension of the system must be 2 to insert a dislocation dipole.")

        if x1==None:
            x1 = self.xmax/3
        if y1==None:
            y1 = self.ymax/2
        if x2==None:
            x2 = 2*self.xmax/3
        if y2==None:
            y2 = self.ymax/2
        if eta==None:
            eta = self.eta0

        sn = self.dislocation_charges[dislocation_type - 1]

        theta1 = self.calc_angle_field_single_vortex([x1,y1])
        theta2 = self.calc_angle_field_single_vortex([x2,y2])

        for n in range(self.number_of_reciprocal_modes):
            if sn[n] != 0:
                eta[n] *= np.exp(1j*sn[n]*theta1)
                eta[n] *= np.exp(-1j*sn[n]*theta2)

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

        # TODO: Implement the PFCper class

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
        self.a = a0*np.array([[1, 0], [1 / 2, np.sqrt(3) / 2], [1 / 2, -np.sqrt(3) / 2]])
        self.q = np.array([[np.sqrt(3) / 2, -1 / 2], [0, 1], [-np.sqrt(3) / 2, -1 / 2]])


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
        self.eta0 = [self.A, self.A, self.A]

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

        a0 = 2 * np.pi
        self.a = a0 * np.array([[1, 0], [0, 1]])
        self.q = np.array([[1, 0], [0, 1], [1, -1], [1, 1]])


        # First initialize the BaseSystem class
        super().__init__(dimension)

        # Type of the system
        self.type = 'PFCsq'

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)



class PFCbcc(PFC):
    def __init__(self, dimension, x_resolution, **kwargs):
        """
        Nothing here yet
        """

        a0 = 2 * np.pi * np.sqrt(2)
        self.a = a0 / 2 * np.array([[-1, 1, 1],
                                    [1, -1, 1],
                                    [1, 1, -1],
                                    [1, 1, 1]])

        self.q = np.array([[0, 1, 1],
                           [1, 0, 1],
                           [1, 1, 0],
                           [0, -1, 1],
                           [-1, 0, 1],
                           [-1, 1, 0]]) / np.sqrt(2)


        # First initialize the BaseSystem class
        super().__init__(dimension)

        # Type of the system
        self.type = 'PFCbcc'


        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)


class PFCfcc(PFC):
    def __init__(self, dimension, x_resolution, **kwargs):
        """
        PFC class for FCC
        """

        a0 = 2 * np.pi * np.sqrt(3)
        self.a = a0 / 2 * np.array([[0, 1, 1],
                                    [1, 0, 1],
                                    [1, 1, 0],
                                    [0, -1, 1],
                                    [-1, 0, 1],
                                    [-1, 1, 0]])

        self.q = np.array([[-1, 1, 1],
                           [1, -1, 1],
                           [1, 1, -1],
                           [1, 1, 1],
                           [2, 0, 0],
                           [0, 2, 0],
                           [0, 0, 2]]) / np.sqrt(3)


        # First initialize the BaseSystem class
        super().__init__(dimension)

        # Type of the system
        self.type = 'PFCfcc'


        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)


class PFCsc(PFC):
    def __init__(self, dimension, x_resolution, **kwargs):
        """
        Nothing here yet
        """

        a0 = 2 * np.pi
        self.a = a0 * np.array([[1, 0, 0],
                                [0, 1, 0],
                                [0, 0, 1]])

        self.q = np.array([[1, 0, 0],
                           [0, 1, 0],
                           [0, 0, 1],
                           [0, 1, 1],
                           [1, 0, 1],
                           [1, 1, 0],
                           [0, -1, 1],
                           [-1, 0, 1],
                           [-1, 1, 0],
                           [-1, 1, 1],
                           [1, -1, 1],
                           [1, 1, -1],
                           [1, 1, 1]])


        # First initialize the BaseSystem class
        super().__init__(dimension)

        # Type of the system
        self.type = 'PFCsc'



        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)


