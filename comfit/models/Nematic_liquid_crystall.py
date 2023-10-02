import numpy as np
from comfit.core.base_system import BaseSystem

class nematic(BaseSystem):

    def __init__(self, dimension, **kwargs):
        """
        Initializes a system to simulate a (active) nematic liquid crystal

        Parameters:
        - dimension : int
            The dimension of the system.
        - x_resolution : int
            The resolution along the x-axis.
        - kwargs : dict, optional
            Optional keyword arguments to set additional parameters.

        Returns:
        - nematic object
            The system object representing the nematic simulation.

        Example:
        nematic = nematic(3, 100, alpha=0.5)
        Creates a BEC system with 3 dimensions and an x-resolution of 100. The activity alpha is set to 0.5.
        """
        super().__init__(dimension, **kwargs)

        # Type of the system
        self.Q = None
        self.Q_f = None
        self.u = None
        self.u_f = None
        self.type = 'Nematic'

        #defoult parameters
        self.alpha = -1  #defult is an extensile system
        self.K = 1
        self.A = 1
        self.B = 1
        self.Lambda = 0 #flow allignment, not sure if this will be implemented
        self.gamma = 1  # rotational diffusion
        self.Gamma = 0 # friction, note in 3 dim this has to be zero
        self.eta = 1 # viscosity


        for key, value in kwargs.items():
            setattr(self, key, value)

    # Initial condition
    def set_initial_condition_disordered(self, noise_strength=0.01):
        """
        Note noise is here only in angle
        :param noise_strength:
        :return:
        """
        if self.dim == 2:
            self.S0 = np.sqrt(self.B)
            theta_rand = noise_strength*np.random.randn(self.xRes,self.yRes)
            self.Q = np.zeros((self.dim,self.dim,self.xRes,self.yRes))
            self.Q[0][0] = self.S0/2 *np.cos(2*theta_rand)
            self.Q[1][1] = -self.S0/2 *np.cos(2*theta_rand)
            self.Q[1][0] =  self.S0/2 *np.sin(2*theta_rand)
            self.Q[0][1] = self.S0/2 *np.sin(2*theta_rand)

            self.Q_f = np.zeros((self.dim,self.dim,self.xRes,self.yRes))
            self.Q_f[0][0] = np.fft.fftn(self.Q[0][0])
            self.Q_f[1][0] = np.fft.fftn(self.Q[1][0])
            self.Q_f[0][1] = np.fft.fftn(self.Q[0][1])
            self.Q_f[1][1] = np.fft.fftn(self.Q[1][1])
        else:
            raise Exception("not included at the moment")