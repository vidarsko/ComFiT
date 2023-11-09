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
            self.S0 = np.sqrt(self.B)#*noise_strength_S
            theta_rand = noise_strength*np.random.randn(self.xRes,self.yRes)
            self.Q = np.zeros((self.dim,self.dim,self.xRes,self.yRes))
            self.Q[0][0] = self.S0/2 *np.cos(2*theta_rand)
            self.Q[1][1] = -self.S0/2 *np.cos(2*theta_rand)
            self.Q[1][0] =  self.S0/2 *np.sin(2*theta_rand)
            self.Q[0][1] = self.S0/2 *np.sin(2*theta_rand)

            self.Q_f = np.fft.fft2(self.Q)

            self.k2 = self.calc_k2() # k2
            self.k2_press = self.calc_k2()
            self.k2_press[0,0] = 1 # for calculating pressure and velocity

        else:
            raise Exception("not included at the moment")

    def calc_evolution_integrating_factors_nematic_f(self):
        """

        :return:  integrating_factors_f
        """


        omega_f = (self.A*self.B-self.K*self.k2  )/self.gamma

        integrating_factors_f = [0, 0, 0]

        integrating_factors_f[0] = np.exp(omega_f * self.dt)
        If1 = integrating_factors_f[0]

        integrating_factors_f[1] = (If1 - 1) / omega_f

        integrating_factors_f[2] = 1 / (self.dt * omega_f ** 2) * (If1 - 1 - omega_f * self.dt)

        return integrating_factors_f

    def calc_nonlinear_evolution_term_no_flow_f(self,Q):
        Q2 = np.sum(Q[i][j]*Q[j][i] for j in range(self.dim) for i in range(self.dim))

        return -2*self.A*np.fft.fftn(Q2 *Q,axes =(range(-self.dim,0)))/self.gamma

    def evolve_nematic_no_flow(self,number_of_steps):
        omega_f = (self.A * self.B - self.K * self.k2) / self.gamma
        integrating_factors_f = self.calc_evolution_integrating_factors_ETDRK4(omega_f)

        for n in range(number_of_steps):
            self.Q, self.Q_f = self.evolve_ETDRK4_loop(integrating_factors_f,self.calc_nonlinear_evolution_term_no_flow_f,
                                                           self.Q, self.Q_f)
        self.Q = np.real(self.Q)

    def calc_S(self):
        if self.dim == 2:
            return 2*np.sqrt((self.Q[0][0])**2 +(self.Q[0][1])**2)

    def calc_u(self,Q):
        '''
        calculate the velocity and its fourier transform. Because of the possibility Gamma = 0 we have to use the self.k2_press to avoid
        division by zero. This is not a problem since the zero mode of all the forces are zero
        :return:
        '''
        self.F_af = self.calc_activ_force_f(Q)
        self.F_pf = self.calc_passive_force_f(Q)
        self.p_f = self.calc_pressure_f()
        grad_pf = self.calc_grad_p_f()
        self.u_f = (self.F_af + self.F_pf-grad_pf )/ (self.Gamma +self.eta*self.k2_press)

        self.u = np.real(np.fft.ifftn(self.u_f, axes=(range(-self.dim, 0))))

    def calc_activ_force_f(self,Q):
        '''
        Function that calculates the activ force in Fourier space.
        :return: active force
        '''
        F_af = []
        for j in range(self.dim):
            F_af.append(np.sum(1j*self.k[i]*np.fft.fftn(self.alpha *Q[j][i]) for i in range(self.dim)))
        return np.array(F_af)

    def calc_passive_force_f(self,Q):
        '''
        Calculate the passive force in Fourier space
        :return: passive force
        '''
        Pi_f = self.calc_passive_stress_f(Q)
        F_pf = []
        for j in range(self.dim):
            F_pf.append(np.sum(1j * self.k[i] *Pi_f[j][i] for i in range(self.dim)))
        return F_pf

    def calc_passive_stress_f(self,Q):
        H = self.calc_molecular_field(Q)
        Antisym_QH = np.zeros_like(self.Q_f)
        Ericksen = np.zeros_like(self.Q_f)
        for i in range(self.dim):
            for j in range(self.dim):
                Antisym_QH[i][j] = np.sum(Q[i][k]*H[k][j] -H[i][k]*Q[k][j] for k in range(self.dim))
                Ericksen[i][j] = - self.K*np.sum(np.fft.ifftn(1j*self.k[i]*np.fft.fftn(Q[m][l]))*
                                          np.fft.ifftn(1j * self.k[j] * np.fft.fftn(Q[m][l]))
                                          for m in range(self.dim) for l in range(self.dim))
        return np.fft.fftn(Ericksen +Antisym_QH, axes=(range(-self.dim, 0)) )
    def calc_molecular_field(self,Q):
        Q2 =  np.sum(Q[i][j]*Q[j][i] for j in range(self.dim) for i in range(self.dim))
        temp = -self.K * np.fft.ifftn( self.k2* np.fft.fftn(Q,axes=(range(-self.dim,0))),axes=(range(-self.dim,0)) )
        return temp +self.A*self.B*Q -2*self.A*Q2*Q
    def calc_pressure_f(self):
        '''
        calculates the pressure in Fourier space. The zero mode is set to zero
        :return: pressure
        '''
        p_af = np.sum(1j*self.k[i]*self.F_af[i] for i in range(self.dim))
        p_pf = np.sum(1j*self.k[i]*self.F_pf[i] for i in range(self.dim))
        return -(p_af + p_pf)/self.k2_press

    def calc_grad_p_f(self):
        grad_pf = []
        for i in range(self.dim):
            grad_pf.append(1j*self.k[i]*self.p_f)
        return np.array(grad_pf)

    def calc_vorticity_tensor(self):
        Omega_f = np.zeros_like(self.Q_f)
        for i in range(self.dim):
            for j in range(self.dim):
                Omega_f[i][j] = (1j*self.k[i]*self.u_f[j] -1j*self.k[j]*self.u_f[i])/2
        Omega = np.fft.ifftn(Omega_f,axes=range(-self.dim,0))
        return np.real(Omega)

    def calc_strain_rate_tensor_f(self):
        E_f = np.zeros_like(self.Q_f)
        for i in range(self.dim):
            for j in range(self.dim):
                E_f[i][j]= (1j*self.k[i]*self.u_f[j] +1j*self.k[j]*self.u_f[i])/2
        return E_f

    def calc_nonlinear_evolution_term_f(self,Q):
        # TODO test and make sure that the passive stress works as intended
        self.calc_u(Q)
        Q_f = np.fft.fftn(Q,axes=range(-self.dim,0))
        N_f = self.calc_nonlinear_evolution_term_no_flow_f(Q)
        Omega =self.calc_vorticity_tensor()
        Antisym_Omega_Q = np.zeros_like(Q_f)

        for i in range(self.dim):
            for j in range(self.dim):
                Antisym_Omega_Q[i][j] = np.sum(Q[i][k]*Omega[k][j] -Omega[i][k]*Q[k][j] for k in range(self.dim))
        advectiv_deriv = - np.sum(self.u[k]* np.fft.ifftn(1j*self.k[k] * Q_f,axes=(range(-self.dim,0)))for k in range(self.dim) )
        return np.fft.fftn(Antisym_Omega_Q +advectiv_deriv, axes=range(-self.dim,0)) +N_f

    def evolve_nematic(self, number_of_steps):
        omega_f = (self.A * self.B - self.K * self.k2) / self.gamma
        integrating_factors_f = self.calc_evolution_integrating_factors_ETDRK4(omega_f)

        for n in range(number_of_steps):
            self.Q, self.Q_f = self.evolve_ETDRK4_loop(integrating_factors_f,
                                                       self.calc_nonlinear_evolution_term_f,
                                                       self.Q, self.Q_f)
        self.Q = np.real(self.Q)


    def calc_defect_density_nematic(self):
        psi0 = np.sqrt(self.B)/2
        psi =[self.Q[0][0],self.Q[0][1]]
        return self.calc_defect_density(psi,psi0)

    def calc_director(self):
        if self.dim == 2:
            psi_n = self.Q[0][0] + 1j*self.Q[0][1]
            angle = np.angle(psi_n)
            return [np.cos(angle/2),np.sin(angle/2)]