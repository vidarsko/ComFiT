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
        Creates a nematic system with 3 dimensions and an x-resolution of 100. The activity alpha is set to 0.5.
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
    def set_initial_condition_disordered(self, noise_strength=0.01,noise_strength_S =0.1):
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

            self.Q_f = np.zeros((self.dim,self.dim,self.xRes,self.yRes),dtype=np.cdouble)
            self.Q_f[0][0] = np.fft.fft2(self.Q[0][0])
            self.Q_f[1][0] = np.fft.fft2(self.Q[1][0])
            self.Q_f[0][1] = np.fft.fft2(self.Q[0][1])
            self.Q_f[1][1] = np.fft.fft2(self.Q[1][1])
        else:
            raise Exception("not included at the moment")

    def calc_evolution_integrating_factors_nematic_f(self):
        """

        :return:  integrating_factors_f
        """
        k2 = self.calc_k2()
        print(self.K)
        omega_f = (self.A*self.B-self.K*k2  )/self.gamma

        integrating_factors_f = [0, 0, 0]

        integrating_factors_f[0] = np.exp(omega_f * self.dt)
        If1 = integrating_factors_f[0]

        integrating_factors_f[1] = (If1 - 1) / omega_f

        integrating_factors_f[2] = 1 / (self.dt * omega_f ** 2) * (If1 - 1 - omega_f * self.dt)

        return integrating_factors_f

    def calc_nonlinear_evolution_term_no_flow_f(self,Q):
        Q2 = np.sum(self.Q[i][j]*self.Q[j][i] for j in range(self.dim) for i in range(self.dim))

        return -2*self.A*np.fft.fftn(Q2 *self.Q,axes =(range(-self.dim,0)))

    def evolve_nematic_no_flow(self,number_of_steps):

        integrating_factors_f = self.calc_evolution_integrating_factors_nematic_f()

        for n in range(number_of_steps):
            self.Q, self.Q_f = self.evolve_ETDRK2_loop(integrating_factors_f,self.calc_nonlinear_evolution_term_no_flow_f,
                                                           self.Q, self.Q_f)
            self.Q = np.real(self.Q)
    def calc_S(self):
        if self.dim == 2:
            return 2*np.sqrt((self.Q[0][0])**2 +(self.Q[0][1])**2)