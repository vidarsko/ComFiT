import numpy
import numpy as np
from comfit.core.base_system import BaseSystem
import scipy as sp
from tqdm import tqdm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import matplotlib.cm as cm
from matplotlib.colors import Normalize

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
        self.Gamma = 0 if 'Gamma' not in kwargs else kwargs['Gamma'] # friction,
        self.eta = 1 if 'eta' not in kwargs else kwargs['eta'] # viscosity


        for key, value in kwargs.items():
            setattr(self, key, value)

    #TODO: Generalize to 3D
    ### defining get functions
    def get_sym(self,Q,i,j):
        """
        Function to axes the tensor element Q_ij of a trace less symmetric matrix, which we will
        save as a vector
        Args:
            i (int): row index
            j (int): colum index
        Returns:
            (numpy.ndarray) Q_ij
        """
        if self.dim == 2:
            return (-1)**(i*j) *Q[(i+j)%2]

        elif self.dim ==3:
            if i != j:
                return Q[i*j + max(i,j)]
            else:
                if i !=2:
                    return Q[3*i]
                else:
                    return -1 * (Q[0] + Q[3])
### TODO see if this can be improved

    def get_anti_sym(self,omega,i,j):
        """
        Function to axes the tensor element omega_ij of an anti-symetric matrix, which we will
                save as a scalar field in 2D and a vector field in 3D
                Args:
                    i (int): row index
                    j (int): colum index
                Returns:
                    (numpy.ndarray) Q_ij
        """
        if self.dim == 2:
            if i==j:
                return 0
            return (-1)**i *omega
        elif self.dim ==3:
            if i ==j:
                return 0
            else:
                return np.sign(j-i)*omega[i+j-1]


    # Initial condition
    def conf_initial_condition_disordered(self, noise_strength=0.01):
        """
        Initialises the system with the nematogens pointing in the x-direction in 2D and in the z-direction in 3D
        with some random noise in the angle.
        Args:
             noise_strength (float): A meshure for how much noise to put in the angle
        returns:
            Initialises self.Q  and self.Q_f
        raises:
            Exception if the dimension is not 2 or 3
        """
        if self.dim == 2:
            self.S0 = np.sqrt(self.B)
            theta_rand = noise_strength*np.random.randn(self.xRes,self.yRes)
            self.Q = np.zeros((2,self.xRes,self.yRes))
            self.Q[0] = self.S0/2 *np.cos(2*theta_rand)
            self.Q[1] =  self.S0/2 *np.sin(2*theta_rand)


            self.Q_f = sp.fft.fft2(self.Q)

            self.k2 = self.calc_k2() # k2
            self.k2_press = self.calc_k2()
            self.k2_press[0,0] = 1 # for calculating pressure and velocity

        elif self.dim == 3:

            #TODO find equilibrium S in three dim for C != 0 (11 /01 /24)

            self.S0 = np.sqrt(3*self.B) / 2

            theta_rand = noise_strength*np.random.randn(self.xRes,self.yRes,self.zRes)
            phi_rand = noise_strength*np.random.randn(self.xRes,self.yRes,self.zRes)

            nx = np.cos(theta_rand)*np.sin(phi_rand)
            ny = np.sin(theta_rand)*np.sin(phi_rand)
            nz = np.cos(theta_rand)

            self.Q = np.zeros((5, self.xRes, self.yRes, self.zRes))
            self.Q[0] = self.S0 *(nx*nx -1/3)
            self.Q[1] = self.S0 *(nx*ny)
            self.Q[2] = self.S0 *(nx*nz)
            self.Q[3] = self.S0 *(ny*ny -1/3)
            self.Q[4] = self.S0 *(ny*nz)

            self.Q_f = sp.fft.fftn(self.Q, axes=(range(-self.dim, 0)))

            self.k2 = self.calc_k2()  # k2
            self.k2_press = self.calc_k2()
            self.k2_press[0, 0, 0] = 1  # for calculating pressure and velocity

        else:
            raise Exception("This dimension is not included for the moment")

    def conf_insert_disclination_dipole(self, dipole_vector=None, dipole_position=None):
        """
        Sets the initial condition for a disclination dipole configuration in a 2-dimensional system.

        Parameters:
            None

        Returns:
            None

        Raises:
            Exception: If the dimension of the system is not 2.
        """
        if not (self.dim == 2):
            raise Exception("The dimension of the system must be 2 for a vortex dipole configuration.")

        if self.Q is None:
            self.conf_initial_condition_disordered(noise_strength=0)

        if dipole_vector is None:
            dipole_vector = [self.xmax / 3, 0]

        if dipole_position is None:
            dipole_position = self.rmid

        psi = self.Q[0] * np.exp(1j * self.calc_angle_field_vortex_dipole(dipole_vector, dipole_position))
        self.Q[0] = np.real(psi)
        self.Q[1] = np.imag(psi)
        self.Q_f = sp.fft.fft2(self.Q)



    def conf_active_channel(self,width,d=7):
        """
        Set the activity to zero everywhere exept for inside a channel of width "width"
        Args:
            width (float): width of the channel
            d (float, optional): width of interface
        returns:
            updates the activity to the channel configuration.
        """
        if self.dim ==2:
            X, Y = np.meshgrid(self.x, self.y, indexing='ij')
            alpha_0 = self.alpha
            self.alpha = alpha_0*(1- 1 / 2 * (2 + np.tanh((X - self.xmid - width/2) / d) - np.tanh((X - self.xmid + width/2) / d)))
        else:
            raise Exception("The active channel is only permitted in three dimensions")

#### calculations related to the flow field
    def conf_u(self,Q):
        #TODO: This function needs a more descriptive name (Vidar 03.01.24)
        '''
        calculate the velocity and its fourier transform. Because of the possibility Gamma = 0 we have to use the self.k2_press to avoid
        division by zero. This is not a problem since the zero mode of all the forces are zero
        :return:
            (numpy.narray) velocity
        '''
        F_af = self.calc_active_force_f(Q)
        F_pf = self.calc_passive_force_f(Q)
        p_f = self.calc_pressure_f(F_af,F_pf)
        grad_pf = self.calc_grad_p_f(p_f)
        self.u_f = (F_af + F_pf-grad_pf )/ (self.Gamma +self.eta*self.k2_press)

        self.u = np.real(sp.fft.ifftn(self.u_f, axes=(range(-self.dim, 0))))

    def calc_active_force_f(self,Q):
        '''
        Function that calculates the activ force in Fourier space.
        :return: active force
        '''
        F_af = []
        for j in range(self.dim):
            F_af.append(np.sum(1j*self.k[i]*sp.fft.fftn(self.alpha *self.get_sym(Q,j,i),axes=(range(-self.dim, 0)) ) for i in range(self.dim)))
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
        return numpy.array(F_pf)

    def calc_passive_stress_f(self,Q):
        #TODO Make 3d
        """
        Calculates the passive stress in fourier space
        Args:
            Q (numpy.narray) the order parameter that we use to find the stress.
        return: (numpy.narray) the passive stress in fourier space
        """
        if self.dim == 2:
            H = self.calc_molecular_field(Q)
            Antisym_QH = np.sum(self.get_sym(Q,0,k)*self.get_sym(H,k,1) -self.get_sym(H,0,k)*self.get_sym(Q,k,1) for k in range(self.dim))
            Ericksen = np.zeros((self.dim,self.xRes,self.yRes),dtype=np.complex128)
            Ericksen[0] = - self.K*np.sum(sp.fft.ifftn(1j*self.k[0]*sp.fft.fftn(self.get_sym(Q,m,l)))*
                                              sp.fft.ifftn(1j * self.k[0] * sp.fft.fftn(self.get_sym(Q,m,l)))
                                              for m in range(self.dim) for l in range(self.dim))
            Ericksen[1] = - self.K*np.sum(sp.fft.ifftn(1j*self.k[0]*sp.fft.fftn(self.get_sym(Q,m,l)))*
                                              sp.fft.ifftn(1j * self.k[1] * sp.fft.fftn(self.get_sym(Q,m,l)))
                                              for m in range(self.dim) for l in range(self.dim))

            stress = np.zeros((self.dim,self.dim,self.xRes,self.yRes),dtype=np.complex128)
            for i in range(self.dim):
                for j in range(self.dim):
                    stress[i][j] = self.get_sym(Ericksen,i,j) + self.get_anti_sym(Antisym_QH,i,j)
            return sp.fft.fftn(stress, axes=(range(-self.dim, 0)) )

        elif self.dim == 3:
            H = self.calc_molecular_field(Q)

            Antisym_QH = np.zeros((3, self.xRes, self.yRes,self.zRes), dtype=np.complex128)
            Antisym_QH[0] = np.sum(self.get_sym(Q,0,k)*self.get_sym(H,k,1) -self.get_sym(H,0,k)*self.get_sym(Q,k,1) for k in range(self.dim))
            Antisym_QH[1] = np.sum(
                self.get_sym(Q, 0, k) * self.get_sym(H, k, 2) - self.get_sym(H, 0, k) * self.get_sym(Q, k, 2) for k in
                range(self.dim))
            Antisym_QH[2] = np.sum(
                self.get_sym(Q, 1, k) * self.get_sym(H, k, 2) - self.get_sym(H, 1, k) * self.get_sym(Q, k, 2) for k in
                range(self.dim))

            Ericksen = np.zeros((5, self.xRes, self.yRes,self.zRes), dtype=np.complex128)
            Ericksen[0] = - self.K * np.sum(sp.fft.ifftn(1j * self.k[0] * sp.fft.fftn(self.get_sym(Q, m, l))) *
                                            sp.fft.ifftn(1j * self.k[0] * sp.fft.fftn(self.get_sym(Q, m, l)))
                                            for m in range(self.dim) for l in range(self.dim))
            Ericksen[1] = - self.K * np.sum(sp.fft.ifftn(1j * self.k[0] * sp.fft.fftn(self.get_sym(Q, m, l))) *
                                            sp.fft.ifftn(1j * self.k[1] * sp.fft.fftn(self.get_sym(Q, m, l)))
                                            for m in range(self.dim) for l in range(self.dim))
            Ericksen[2] = - self.K * np.sum(sp.fft.ifftn(1j * self.k[0] * sp.fft.fftn(self.get_sym(Q, m, l))) *
                                            sp.fft.ifftn(1j * self.k[2] * sp.fft.fftn(self.get_sym(Q, m, l)))
                                            for m in range(self.dim) for l in range(self.dim))
            Ericksen[3] = - self.K * np.sum(sp.fft.ifftn(1j * self.k[1] * sp.fft.fftn(self.get_sym(Q, m, l))) *
                                            sp.fft.ifftn(1j * self.k[1] * sp.fft.fftn(self.get_sym(Q, m, l)))
                                            for m in range(self.dim) for l in range(self.dim))
            Ericksen[4] = - self.K * np.sum(sp.fft.ifftn(1j * self.k[1] * sp.fft.fftn(self.get_sym(Q, m, l))) *
                                            sp.fft.ifftn(1j * self.k[2] * sp.fft.fftn(self.get_sym(Q, m, l)))
                                            for m in range(self.dim) for l in range(self.dim))

            stress = np.zeros((self.dim,self.dim,self.xRes,self.yRes,self.zRes))
            for i in range(self.dim):
                for j in range(self.dim):
                    stress[i][j] = self.get_sym(Ericksen,i,j) + self.get_anti_sym(Antisym_QH,i,j)
            return sp.fft.fftn(stress, axes=(range(-self.dim, 0)) )


    def calc_molecular_field(self,Q):
        """
        Finds the molecular field (NB! need to be rewriten when C != 0)
        Args
            Q (numpy.ndarray): The nematic tensor
        Returns:
             (numpy.ndarray): The molecular field
        """
        #TODO make this compatible with C != 0 (3D ok when C==0)
        Q2 =  np.sum(self.get_sym(Q,i,j)*self.get_sym(Q,j,i) for j in range(self.dim) for i in range(self.dim))
        temp = -self.K * sp.fft.ifftn( self.k2* sp.fft.fftn(Q,axes=(range(-self.dim,0))),axes=(range(-self.dim,0)) )
        return temp +self.A*self.B*Q -2*self.A*Q2*Q

    def calc_pressure_f(self,F_af,F_pf):
        '''
        calculates the pressure in Fourier space. The zero mode is set to zero
        Returns:
            (numpy.ndarray) the pressure
        '''
        p_af = np.sum(1j*self.k[i]*F_af[i] for i in range(self.dim))
        p_pf = np.sum(1j*self.k[i]*F_pf[i] for i in range(self.dim))
        return -(p_af + p_pf)/self.k2_press

    def calc_grad_p_f(self,p_f):
        """
        Caclulates the gradient of the pressure
        Returns:
             (numpy.ndarray) gradient of th epressure
        """
        grad_pf = []
        for i in range(self.dim):
            grad_pf.append(1j*self.k[i]*p_f)
        return np.array(grad_pf)

    def calc_vorticity_tensor(self):
        """
        Calculates the vorticity tensor
        returns:
            (numpy.ndarray) The vorticity tensor
        """
        if self.dim == 2:

            Omega_f = (1j*self.k[0]*self.u_f[1] -1j*self.k[1]*self.u_f[0])/2
            Omega = sp.fft.ifftn(Omega_f,axes=range(-self.dim,0))
            return np.real(Omega)

        elif self.dim ==3:
            Omega = np.zeros((3,self.xRes,self.yRes,self.zRes))
            Omega[0] = np.real(sp.fft.ifftn((1j*self.k[0]*self.u_f[1] -1j*self.k[1]*self.u_f[0])/2,axes=range(-self.dim,0)))
            Omega[1] = np.real(sp.fft.ifftn((1j * self.k[0] * self.u_f[2] - 1j * self.k[2] * self.u_f[0]) / 2,
                                            axes=range(-self.dim, 0)))
            Omega[2] = np.real(sp.fft.ifftn((1j * self.k[1] * self.u_f[2] - 1j * self.k[2] * self.u_f[1]) / 2,
                                            axes=range(-self.dim, 0)))
            return Omega

    def calc_strain_rate_tensor_f(self):
        # TODO Make symetric and 3d
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
        #TODO make 3d compatible
        """
        Calculates the non-linear evolution function for the nematic
        Args:
            Q (numpy.narray) the nematic orderparameter
        returns:
            (numpy.narray) the non-linear evolution function evaluated in Fourier space
        """
        if self.dim == 2:
            self.conf_u(Q)
            Q_f = sp.fft.fftn(Q,axes=range(-self.dim,0))
            N_f = self.calc_nonlinear_evolution_term_no_flow_f(Q)
            Omega =self.calc_vorticity_tensor()
            Antisym_Omega_Q = np.zeros_like(Q_f)

            Antisym_Omega_Q[0] = np.sum(self.get_sym(Q,0,k)*self.get_anti_sym(Omega,k,0) -
                                        self.get_anti_sym(Omega,0,k)*self.get_sym(Q,k,0) for k in range(self.dim))
            Antisym_Omega_Q[1] = np.sum(
                self.get_sym(Q, 0, k) * self.get_anti_sym(Omega,k,1) - self.get_anti_sym(Omega,0,k) * self.get_sym(Q, k, 1) for k in range(self.dim))
            advectiv_deriv = - np.sum(self.u[k]* sp.fft.ifftn(1j*self.k[k] * Q_f,axes=(range(-self.dim,0)))for k in range(self.dim) )
            return sp.fft.fftn(Antisym_Omega_Q +advectiv_deriv, axes=range(-self.dim,0)) +N_f

        elif self.dim == 3:
            self.conf_u(Q)
            Q_f = sp.fft.fftn(Q, axes=range(-self.dim, 0))
            N_f = self.calc_nonlinear_evolution_term_no_flow_f(Q)
            Omega = self.calc_vorticity_tensor()

            advectiv_deriv = - np.sum(
                self.u[k] * sp.fft.ifftn(1j * self.k[k] * Q_f, axes=(range(-self.dim, 0))) for k in range(self.dim))

            Antisym_Omega_Q = np.zeros_like(Q_f)

            Antisym_Omega_Q[0] = np.sum(self.get_sym(Q, 0, k) * self.get_anti_sym(Omega, k, 0) -
                                        self.get_anti_sym(Omega, 0, k) * self.get_sym(Q, k, 0) for k in range(self.dim))
            Antisym_Omega_Q[1] = np.sum(
                self.get_sym(Q, 0, k) * self.get_anti_sym(Omega, k, 1) - self.get_anti_sym(Omega, 0, k) * self.get_sym(
                    Q, k, 1) for k in range(self.dim))
            Antisym_Omega_Q[2] = np.sum(
                self.get_sym(Q, 0, k) * self.get_anti_sym(Omega, k, 2) - self.get_anti_sym(Omega, 0, k) * self.get_sym(
                    Q, k, 2) for k in range(self.dim))
            Antisym_Omega_Q[3] = np.sum(
                self.get_sym(Q, 1, k) * self.get_anti_sym(Omega, k, 1) - self.get_anti_sym(Omega, 1, k) * self.get_sym(
                    Q, k, 1) for k in range(self.dim))
            Antisym_Omega_Q[4] = np.sum(
                self.get_sym(Q, 1, k) * self.get_anti_sym(Omega, k, 2) - self.get_anti_sym(Omega, 1, k) * self.get_sym(
                    Q, k, 2) for k in range(self.dim))
            return sp.fft.fftn(Antisym_Omega_Q + advectiv_deriv, axes=range(-self.dim, 0)) + N_f

        else:
            raise Exception("This dimension is not implemented at the moment")

    def calc_nonlinear_evolution_term_no_flow_f(self,Q):
        #TODO generalize to 3d (think it alredy is for C==0)
        """
            Calculates the non-linear evolution function for the nematic without the flow field
                Args:
                    Q (numpy.narray) the nematic orderparameter
                returns:
                    (numpy.narray) the non-linear evolution function evaluated in Fourier space
                """
        Q2 = np.sum(self.get_sym(Q,i,j)*self.get_sym(Q,j,i) for j in range(self.dim) for i in range(self.dim))

        return -2*self.A*sp.fft.fftn(Q2 *Q,axes =(range(-self.dim,0)))/self.gamma


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
        for n in tqdm(range(number_of_steps), desc='evolving the active nematic'):
            self.Q, self.Q_f = solver(integrating_factors_f,
                                                       self.calc_nonlinear_evolution_function_f,
                                                       self.Q, self.Q_f)
        self.Q = np.real(self.Q)

    def evolve_nematic_no_flow(self,number_of_steps,method = 'ETD2RK'):
        ''' Evolver for the nematic system without the flow field
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

        for n in tqdm(range(number_of_steps), desc='evolving the active nematic'):
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
        if self.dim == 2:
            psi0 = np.sqrt(self.B)/2
            psi =[self.Q[0],self.Q[1]]
        else:
            raise Exception("Not implemented")
        return self.calc_defect_density(psi,psi0)

    def calc_dt_psi(self,Q_prev,delta_t):
        dt_Q = (self.Q -Q_prev)/delta_t
        return dt_Q[0] + 1j*dt_Q[1]

    def calc_S(self):
        # TODO: This function will be removed soon
        '''
        Calculates the strength of nematic order S
        :returns
            (numpy.narray) S
        '''
        if self.dim == 2:
            return 2 * np.sqrt((self.Q[0]) ** 2 + (self.Q[1]) ** 2)

        elif self.dim == 3:
            Q2 = np.sum(self.get_sym(self.Q, i, j) * self.get_sym(self.Q, j, i) for j in range(self.dim) for i in
                        range(self.dim))
            return np.sqrt(3 * Q2 / 2)

    def calc_order_and_director(self):
        """
        Finds the amount of order (S) and the director field (n)
        return:
            (tuple): (scalar field) S - amount of order   , (vector field) the director field
        """
        if self.dim == 2:
            psi_n = self.Q[0] + 1j*self.Q[1]
            angle = np.angle(psi_n)
            S =2 * np.sqrt((self.Q[0]) ** 2 + (self.Q[1]) ** 2)
            return S, [np.cos(angle/2),np.sin(angle/2)]
        elif self.dim ==3:
        ## need to construct a matrix
            Q_eig = numpy.zeros((self.xRes,self.yRes,self.zRes,self.dim,self.dim))
            for i in range(self.dim):
                for j in range(self.dim):
                    Q_eig[:,:,:,i,j] = self.get_sym(self.Q,i,j)

            eigvals, eigvectors = numpy.linalg.eigh(Q_eig)
            S = 3/2 *eigvals[:,:,:,2]
            n = np.transpose(eigvectors[:,:,:,2], (3,0,1,2))

            return S, n

    def calc_vortex_velocity_field(self, dt_Q, psi=None):
        """
        Calculates the velocity field of the defect in two dimensions 
        :param dt_psi:
        :param psi:
        :return:
        """
        if self.dim ==2:
            if psi is None:
                psi = self.Q[0] +1j*self.Q[1]

            dt_psi = dt_Q[0] + 1j * dt_Q[1]

            return self.calc_defect_velocity_field([np.real(psi), np.imag(psi)],
                                                   [np.real(dt_psi), np.imag(dt_psi)])

    def calc_defect_polarization_field(self):
        ex = np.real(sp.fft.ifftn(1j*self.k[0]*self.Q_f[0] + 1j*self.k[1]*self.Q_f[1]))
        ey = np.real(sp.fft.ifftn(1j * self.k[0] * self.Q_f[1] - 1j * self.k[1] * self.Q_f[0]))
        return np.array([ex,ey])

    def calc_disclination_nodes_nem(self, dt_Q=None,polarization = None):
        """
        Calculate the positions and charges of vortex nodes based on the defect density.
        Returns:
            list: A list of dictionaries representing the vortex nodes. Each dictionary contains the following keys:
                  - 'position_index': The position index of the vortex node in the defect density array.
                  - 'charge': The charge of the vortex node.
                  - 'position': The position of the vortex node as a list [x, y].
                  - 'velocity': The velocity of the vortex node as a list [vx, vy].
        """

        # Calculate defect density
        if self.dim == 2:
            psi = self.Q[0] + 1j*self.Q[1]
            rho = self.calc_defect_density_nematic()

            if dt_Q is not None:
                velocity_field = self.calc_vortex_velocity_field(dt_Q, psi)


            vortex_nodes = self.calc_defect_nodes(np.abs(rho))
            for vortex in vortex_nodes:
                vortex['charge'] = np.sign(rho[vortex['position_index']])

                if dt_Q is not None:
                    vortex['velocity'] = [velocity_field[0][vortex['position_index']],
                                          velocity_field[1][vortex['position_index']]]
                else:
                    vortex['velocity'] = [float('nan'), float('nan')]

                if vortex['charge'] >0 and polarization is not None:
                    vortex['polarization'] = [polarization[0][vortex['position_index']]
                        ,polarization[1][vortex['position_index']]]
                else:
                    vortex['polarization'] = [float('nan'), float('nan')]


        return vortex_nodes

    def plot_disclination_nodes(self, vortex_nodes, ax=None):

        if self.dim == 2:

            if ax == None:
                ax = plt.gcf().add_subplot(111)

            x_coords_pos = []
            y_coords_pos = []

            x_coords_neg = []
            y_coords_neg = []

            vx_coords_pos = []
            vy_coords_pos = []


            vx_coords_neg = []
            vy_coords_neg = []

            px_coords_pos = []
            py_coords_pos = []



            for vortex in vortex_nodes:

                if vortex['charge'] > 0:
                    x_coords_pos.append(vortex['position'][0])
                    y_coords_pos.append(vortex['position'][1])
                    vx_coords_pos.append(vortex['velocity'][0])
                    vy_coords_pos.append(vortex['velocity'][1])
                    px_coords_pos.append(vortex['polarization'][0])
                    py_coords_pos.append(vortex['polarization'][1])
                else:
                    x_coords_neg.append(vortex['position'][0])
                    y_coords_neg.append(vortex['position'][1])
                    vx_coords_neg.append(vortex['velocity'][0])
                    vy_coords_neg.append(vortex['velocity'][1])


            # print(x_coords_pos,y_coords_pos)
            # print(x_coords_neg,y_coords_neg)
            ax.scatter(x_coords_pos, y_coords_pos, marker='+', color='red')
            ax.scatter(x_coords_neg, y_coords_neg, marker='*', color='blue')
            ax.quiver(x_coords_pos, y_coords_pos, vx_coords_pos, vy_coords_pos, color='black')
            ax.quiver(x_coords_neg, y_coords_neg, vx_coords_neg, vy_coords_neg, color='black')
            ax.quiver(x_coords_pos, y_coords_pos, px_coords_pos, py_coords_pos, color='red')
            ax.set_aspect('equal')
            ax.set_facecolor('none')

            ax.set_xlabel('$x/a_0$')
            ax.set_ylabel('$y/a_0$')

            ax.set_xlim([0, self.xmax-self.dx])
            ax.set_ylim([0, self.ymax-self.dy])


    def plot_field_velocity_and_director(self, field, velocity, director, ax=None, colorbar=True, colormap='viridis',
                                         cmax=None, cmin=None,
                                         number_of_layers=1, hold=False):
        """
        Plot the field, velocity, and director in the given axes.

        Parameters:
            field (ndarray): The field to be plotted.
            velocity (ndarray): The velocity to be plotted.
            director (ndarray): The director to be plotted.
            ax (Axes, optional): The axes to plot the field, velocity, and director on. If not provided, a new subplot will be created.
            colorbar (bool, optional): Whether to show the colorbar. Default is True.
            colormap (str, optional): The colormap to use for plotting the field. Default is 'viridis'.
            cmax (float, optional): The maximum value for the colorbar. If not provided, the maximum value of the field will be used.
            cmin (float, optional): The minimum value for the colorbar. If not provided, the minimum value of the field will be used.
            number_of_layers (int, optional): The number of layers in the plot. Default is 1.
            hold (bool, optional): Whether to hold the plot. Default is False.

        Returns:
            ax (Axes): The axes with the plotted field, velocity, and director.

        Raises:
            Exception: If the plotting function is not yet configured for dimensions other than 2.

        Note: streamplot assumes xy indexing and not ij. I think it is suficient
        just transpose the matrices before putting them in
        """
        if self.dim == 2:
            if ax == None:
                ax = plt.gcf().add_subplot(111)
            X, Y = np.meshgrid(self.x, self.y, indexing='ij')
            cbar = ax.pcolormesh(X, Y, field, shading='gouraud', cmap=colormap)
            plt.colorbar(cbar)
            ax.streamplot(X.T, Y.T, (velocity[0]).T, (velocity[1]).T, color='w')
            ax.quiver(X, Y, director[0], director[1], headwidth=0, scale=50)
            ax.quiver(X, Y, -director[0], -director[1], headwidth=0, scale=50)
            ax.set_aspect('equal')
            return ax

        else:
            raise Exception("This plotting function not yet configured for other dimension")



    def plot_nematic_3D(self,n,S,ax = None,step=None):
        # TODO find a better way of plotting.
        if self.dim != 3:
            raise Exception("Dimension not allowed")
        if ax == None:
            ax = plt.gcf().add_subplot(111, projection='3d')

        if step is None:
            step = 2



        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

        X_plot = X[::step, ::step, ::step]
        Y_plot = Y[::step, ::step, ::step]
        Z_plot = Z[::step, ::step, ::step]
        U_plot = n[0][::step, ::step, ::step]
        V_plot = n[1][::step, ::step, ::step]
        W_plot = n[2][::step, ::step, ::step]
        c = S[::step,::step,::step]

        c  = (c.ravel() - c.min()) / c.ptp()
        c = np.concatenate((c, np.repeat(c, 2)))

        c = plt.cm.viridis(c)





        ax.quiver(X_plot, Y_plot, Z_plot, U_plot, V_plot, W_plot, arrow_length_ratio=0,pivot= "middle",colors=c)

        ax.set_xlabel('$x/a_0$')
        ax.set_ylabel('$y/a_0$')
        ax.set_zlabel('$z/a_0$')
        ax.set_aspect('equal')
        ax.set_xlim([0, self.xmax - self.dx])
        ax.set_ylim([0, self.ymax - self.dy])
        ax.set_zlim([0, self.zmax - self.dz])