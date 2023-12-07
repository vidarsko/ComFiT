import numpy as np
from comfit.core.base_system import BaseSystem

class NematicLiquidCrystal(BaseSystem):

    def __init__(self, dimension, **kwargs):
        """
        Initializes a system to simulate a (active) nematic liquid crystal

        Parameters:
        - dimension : int
            The dimension of the system. Note that only d=2 is implemented at the moment
        - x_resolution : int
            The resolution along the x-axis.
        - kwargs : dict, optional
            Optional keyword arguments to set additional parameters.

        Returns:
        - nematic object
            The system object representing the nematic simulation.

        Example:
        nematic = nematic(2, 100, alpha=-0.5)
        Creates a nematic system with 2 dimensions and an x-resolution of 100. The activity alpha is set to -0.5.
        """
        super().__init__(dimension, **kwargs)

        # Type of the system
        self.Q = None
        self.Q_f = None
        self.u = None
        self.u_f = None
        self.type = 'NematicLiquidCrystal'

        #defoult parameters
        self.alpha = -1 if 'alpha' not in kwargs else kwargs['alpha']   #defult is an extensile system
        self.K = 1 if 'K' not in kwargs else kwargs['K']
        self.A = 1 if 'A' not in kwargs else kwargs['A']
        self.B = 1 if 'B' not in kwargs else kwargs['B']
        self.Lambda = 0 if 'Lambda' not in kwargs else kwargs['Lambda'] #flow allignment, not sure if this will be implemented
        self.gamma = 1  if 'gamma' not in kwargs else kwargs['gamma']  # rotational diffusion
        self.Gamma = 0 if 'Gamma' not in kwargs else kwargs['Gamma'] # friction, note in 3 dim this has to be zero
        self.eta = 1 if 'eta' not in kwargs else kwargs['eta'] # viscosity


        for key, value in kwargs.items():
            setattr(self, key, value)

    # Initial condition
    def conf_initial_condition_disordered(self, noise_strength=0.01):
        """
        Initialises the system with the nematogens pointing in the x-direction with some random noise in the angle
        Args:
             noise_strength (float): A meshure for how much noise to put in the angle
        returns:
            Initialises self.Q  and self.Q_f
        """
        if self.dim == 2:
            self.S0 = np.sqrt(self.B)
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




    def calc_S(self):
        '''
        Calculates the strength of nematic order S
        :returns
            (numpy.narray) S
        '''
        if self.dim == 2:
            return 2*np.sqrt((self.Q[0][0])**2 +(self.Q[0][1])**2)

#### calculations related to the flow field
    def calc_u(self,Q):
        '''
        calculate the velocity and its fourier transform. Because of the possibility Gamma = 0 we have to use the self.k2_press to avoid
        division by zero. This is not a problem since the zero mode of all the forces are zero
        :return:
            (numpy.narray) velocity
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
        """
        Calculates the passive stress in fourier space
        Args:
            Q (numpy.narray) the order parameter that we use to find the stress.
        return: (numpy.narray) the passive stress in fourier space
        """
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

        #TODO: Double check that Q[1][1] is the same as Q[1,1,:,:] (Vidar 16.11.23)
        # Answer: this is tested in 231121jr_test_of_nematic_indexing.py,
        # where I chech that they have the same shape and that the difference between the
        # absolute values are zero (Jonas 21.11.23)

    def calc_molecular_field(self,Q):
        """
        Finds the molecular field (NB! strictly 2D at the moment)
        Args
            Q (numpy.ndarray): The nematic tensor
        Returns:
             (numpy.ndarray): The molecular field
        """
        Q2 =  np.sum(Q[i][j]*Q[j][i] for j in range(self.dim) for i in range(self.dim))
        temp = -self.K * np.fft.ifftn( self.k2* np.fft.fftn(Q,axes=(range(-self.dim,0))),axes=(range(-self.dim,0)) )
        return temp +self.A*self.B*Q -2*self.A*Q2*Q

    def calc_pressure_f(self):
        '''
        calculates the pressure in Fourier space. The zero mode is set to zero
        Returns:
            (numpy.ndarray) the pressure
        '''
        p_af = np.sum(1j*self.k[i]*self.F_af[i] for i in range(self.dim))
        p_pf = np.sum(1j*self.k[i]*self.F_pf[i] for i in range(self.dim))
        return -(p_af + p_pf)/self.k2_press

    def calc_grad_p_f(self):
        """
        Caclulates the gradient of the pressure
        Returns:
             (numpy.ndarray) gradient of th epressure
        """
        grad_pf = []
        for i in range(self.dim):
            grad_pf.append(1j*self.k[i]*self.p_f)
        return np.array(grad_pf)

    def calc_vorticity_tensor(self):
        """
        Calculates the vorticity tensor
        returns:
            (numpy.ndarray) The vorticity tensor
        """
        Omega_f = np.zeros_like(self.Q_f)
        for i in range(self.dim):
            for j in range(self.dim):
                Omega_f[i][j] = (1j*self.k[i]*self.u_f[j] -1j*self.k[j]*self.u_f[i])/2
        Omega = np.fft.ifftn(Omega_f,axes=range(-self.dim,0))
        return np.real(Omega)

    def calc_strain_rate_tensor_f(self):
        """
        Calculates the strainrate tensor
        returns:
            (numpy.ndarray) The strainrate
        """
        E_f = np.zeros_like(self.Q_f)
        for i in range(self.dim):
            for j in range(self.dim):
                E_f[i][j]= (1j*self.k[i]*self.u_f[j] +1j*self.k[j]*self.u_f[i])/2
        return E_f

#### Calculation of non-linear evolution terms
    def calc_nonlinear_evolution_function_f(self,Q):
        # TODO test and make sure that the passive stress works as intended (Jonas: 2023/11/14)
        """
        Calculates the non-linear evolution function for the nematic
        Args:
            Q (numpy.narray) the nematc orderparameter
        returns:
            (numpy.narray) the non-linear evolution function evaluated in Fourier space
        """
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

    def calc_nonlinear_evolution_term_no_flow_f(self,Q):
        """
            Calculates the non-linear evolution function for the nematic without the flow field
                Args:
                    Q (numpy.narray) the nematc orderparameter
                returns:
                    (numpy.narray) the non-linear evolution function evaluated in Fourier space
                """
        Q2 = np.sum(Q[i][j]*Q[j][i] for j in range(self.dim) for i in range(self.dim))

        return -2*self.A*np.fft.fftn(Q2 *Q,axes =(range(-self.dim,0)))/self.gamma


##### evolvers
    def evolve_nematic(self, number_of_steps, method= 'ETD2RK'):
        '''
         Evolver for the nematic system
            Args:
                number_of_steps (int) the number of time steps that we are evolving the equation
                method (string, optional) the integration method we want to use. ETD2RK is sett as default
            returns:
                Updates the fields self.Q and self.Q_f
        '''
        omega_f = (self.A * self.B - self.K * self.k2) / self.gamma

        if method == 'ETD2RK':
            integrating_factors_f = self.calc_evolution_integrating_factors_ETD2RK(omega_f)
            solver = self.evolve_ETD2RK_loop

        elif method == 'ETD4RK':
            integrating_factors_f = self.calc_evolution_integrating_factors_ETD4RK(omega_f)
            solver = self.evolve_ETD4RK_loop
        else:
            raise Exception('Mehtod is not implemented')
        for n in range(number_of_steps):
            self.Q, self.Q_f = solver(integrating_factors_f,
                                                       self.calc_nonlinear_evolution_function_f,
                                                       self.Q, self.Q_f)
        self.Q = np.real(self.Q)

    def evolve_nematic_no_flow(self,number_of_steps,method = 'ETD2RK'):
        '''
                 Evolver for the nematic system without the flow field
                    Args:
                        number_of_steps (int) the number of time steps that we are evolving the equation
                        method (string, optional) the integration method we want to use. ETD2RK is sett as default
                    returns:
                        Updates the fields self.Q and self.Q_f
                '''
        omega_f = (self.A * self.B - self.K * self.k2) / self.gamma

        if method == 'ETD2RK':
            integrating_factors_f = self.calc_evolution_integrating_factors_ETD2RK(omega_f)
            solver = self.evolve_ETD2RK_loop

        elif method == 'ETD4RK':
            integrating_factors_f = self.calc_evolution_integrating_factors_ETD4RK(omega_f)
            solver = self.evolve_ETD4RK_loop
        else:
            raise Exception('Mehtod is not implemented')

        for n in range(number_of_steps):
            self.Q, self.Q_f = solver(integrating_factors_f,self.calc_nonlinear_evolution_term_no_flow_f,
                                                           self.Q, self.Q_f)
        self.Q = np.real(self.Q)


##### defect tracking
    def calc_defect_density_nematic(self):
        """
        calculates the defect density
        return:
            (numpy.narray) The defect density
        """
        psi0 = np.sqrt(self.B)/2
        psi =[self.Q[0][0],self.Q[0][1]]
        return self.calc_defect_density(psi,psi0)

    def calc_director(self):
        """
        Finds the director field
        return:
            (numpy.narray) the director field
        """
        if self.dim == 2:
            psi_n = self.Q[0][0] + 1j*self.Q[0][1]
            angle = np.angle(psi_n)
            return [np.cos(angle/2),np.sin(angle/2)]