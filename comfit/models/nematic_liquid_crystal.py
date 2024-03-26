import numpy
import numpy as np
from comfit.core.base_system import BaseSystem
import scipy as sp
from tqdm import tqdm
import matplotlib.pyplot as plt
from comfit.tools.tool_colormaps import tool_colormap_angle, tool_colormap_bluewhitered, tool_colormap_sunburst
from comfit.tools.tool_create_orthonormal_triad import tool_create_orthonormal_triad

from comfit.tools.tool_math_functions import levi_civita_symbol
from mpl_toolkits.mplot3d import axes3d
import matplotlib.cm as cm
from matplotlib.colors import Normalize

class NematicLiquidCrystal(BaseSystem):

    def __init__(self, dimension, **kwargs):
        """
        Initializes a system to simulate a (active) nematic liquid crystal

        Input:
        - dimension : int
            The dimension of the system. Note that only d=2 is implemented at the moment
        - x_resolution : int
            The resolution along the x-axis.
        - kwargs : dict, optional
            Optional keyword arguments to set additional parameters.

        Output:
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
        self.C = 0 if 'C' not in kwargs  else kwargs['C']  # note: in two dim C is never called
        self.Lambda = 0 if 'Lambda' not in kwargs else kwargs['Lambda'] #flow allignment, not sure if this will be implemented
        self.gamma = 1  if 'gamma' not in kwargs else kwargs['gamma']  # rotational diffusion
        self.Gamma = 0 if 'Gamma' not in kwargs else kwargs['Gamma'] # friction,
        self.eta = 1 if 'eta' not in kwargs else kwargs['eta'] # viscosity


        for key, value in kwargs.items():
            setattr(self, key, value)

    def __str__(self):
        """
        Output a string representation of the system.
        """
        return "NematicLiquidCrystal"

    # Initial condition
    def conf_initial_condition_ordered(self, noise_strength=0.01):
        """
        Initialises the system with the nematogens pointing in the x-direction in 2D and in the z-direction in 3D
        with some random noise in the angle.
        
        Input:
             noise_strength (float): A meshure for how much noise to put in the angle
        
        Output:
            Initialises self.Q  and self.Q_f
        
        Raises:
            Exception if the dimension is not 2 or 3
        """
        if self.dim == 2:
            S0 = np.sqrt(self.B)
            theta_rand = noise_strength*np.random.randn(self.xRes,self.yRes)
            self.Q = np.zeros((2,self.xRes,self.yRes))
            self.Q[0] = S0/2 *np.cos(2*theta_rand)
            self.Q[1] =  S0/2 *np.sin(2*theta_rand)


            self.Q_f = sp.fft.fft2(self.Q)

            self.k2 = self.calc_k2() # k2
            self.k2_press = self.calc_k2()
            self.k2_press[0,0] = 1 # for calculating pressure and velocity

        elif self.dim == 3:

            S0 = 1/8* self.C/self.A + 1/2 * np.sqrt(self.C**2 /(16*self.A**2) + 3*self.B)


            theta_rand = noise_strength*np.random.randn(self.xRes,self.yRes,self.zRes)
            phi_rand = noise_strength*np.random.randn(self.xRes,self.yRes,self.zRes)

            nx = np.cos(theta_rand)*np.sin(phi_rand)
            ny = np.sin(theta_rand)*np.sin(phi_rand)
            nz = np.cos(phi_rand)

            self.Q = np.zeros((5, self.xRes, self.yRes, self.zRes))
            self.Q[0] = S0 *(nx*nx -1/3)
            self.Q[1] = S0 *(nx*ny)
            self.Q[2] = S0 *(nx*nz)
            self.Q[3] = S0 *(ny*ny -1/3)
            self.Q[4] = S0 *(ny*nz)

            self.Q_f = sp.fft.fftn(self.Q, axes=(range(-self.dim, 0)))

            self.k2 = self.calc_k2()  # k2
            self.k2_press = self.calc_k2()
            self.k2_press[0, 0, 0] = 1  # for calculating pressure and velocity

        else:
            raise Exception("This dimension is not included for the moment")

    def conf_insert_disclination_dipole(self, dipole_vector=None, dipole_position=None):
        """
        Sets the initial condition for a disclination dipole configuration in a 2-dimensional system.

        Input:
            None

        Output:
            None

        Raises:
            Exception: If the dimension of the system is not 2.
        """
        if not (self.dim == 2):
            raise Exception("The dimension of the system must be 2 for a vortex dipole configuration.")

        if self.Q is None:
            self.conf_initial_condition_ordered(noise_strength=0)

        if dipole_vector is None:
            dipole_vector = [self.xmax / 3, 0]

        if dipole_position is None:
            dipole_position = self.rmid

        psi = (self.Q[0] + 1j*self.Q[1]) * np.exp(1j * self.calc_angle_field_vortex_dipole(dipole_vector, dipole_position))
        self.Q[0] = np.real(psi)
        self.Q[1] = np.imag(psi)
        self.Q_f = sp.fft.fft2(self.Q)

    def conf_initial_disclination_lines(self, position=None):
        """
        Sets the initial condition for a disclination line in a 3-dimensional system.
        The dislocation is parralell to the z-axis

        Input:
            position (list): the position of the dislocation. Only the position in the xy plane is used
            angle (float): the angle between the director and the z axis.
            sign (float): +1 or -1. Charge of the dislocation in the xy-plain

        Output:
            Sets the value of self.Q and self.Q_f
        """
        if not (self.dim == 3):
            raise Exception("The dimension of the system must be 3 for a disclination line configuration.")

        if position is None:
            position1 = [self.xmid+self.xmax/3, self.ymid]
            position2 = [self.xmid - self.xmax/3, self.ymid]


        theta_1 = 1/2*np.arctan2((self.y-position1[1]),(self.x-position1[0]))
        theta_2 = -1/2*np.arctan2((self.y-position2[1]),(self.x-position2[0]))
        theta = theta_1 +theta_2

        S0 = 1/8* self.C/self.A + 1/2 * np.sqrt(self.C**2 /(16*self.A**2) + 3*self.B)

        nx =  np.cos(theta)
        ny = np.sin(theta)
        nz = np.zeros_like(nx)

        self.Q = np.zeros((5, self.xRes, self.yRes, self.zRes))
        self.Q[0] = S0 * (nx * nx - 1 / 3)
        self.Q[1] = S0 * (nx * ny)
        self.Q[2] = S0 * (nx * nz)
        self.Q[3] = S0 * (ny * ny - 1 / 3)
        self.Q[4] = S0 * (ny * nz)

        self.Q_f = sp.fft.fftn(self.Q, axes=(range(-self.dim, 0)))

        self.k2 = self.calc_k2()  # k2
        self.k2_press = self.calc_k2()
        self.k2_press[0, 0, 0] = 1




    def conf_active_channel(self,width,d=7):
        """
        Set the activity to zero everywhere exept for inside a channel of width "width"
        
        Input:
            width (float): width of the channel
            d (float, optional): width of interface
        
        Output:
            updates the activity to the channel configuration.
        """
        if self.dim ==2:
            X, Y = np.meshgrid(self.x, self.y, indexing='ij')
            alpha_0 = self.alpha
            self.alpha = alpha_0*(1- 1 / 2 * (2 + np.tanh((X - self.xmid - width/2) / d) - np.tanh((X - self.xmid + width/2) / d)))
        else:
            raise Exception("The active channel is only permitted in two dimensions")

#### calculations related to the flow field
    def conf_u(self,Q):
        #TODO: This function needs a more descriptive name (Vidar 03.01.24)
        '''
        Updates the velocity and its fourier transform given a nematic field Q.
        
        Input:
            (numpy.narray) the Q tensor
        
        Output:
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
        
        Input:
            Q (numpy.narray) the order parameter that we use to find the force.
        
        Output:
            (numpy.narray) the active force in Fourier space
        '''
        F_af = []
        for j in range(self.dim):
            F_af.append(np.sum(1j*self.k[i]*sp.fft.fftn(self.alpha *self.get_sym_tl(Q,j,i),axes=(range(-self.dim, 0)) ) for i in range(self.dim)))
        return np.array(F_af)

    def calc_passive_force_f(self,Q):
        '''
        Calculate the passive force in Fourier space
        
        Input:
            Q (numpy.narray) the order parameter that we use to find the force.
        
        Output: 
            numpy.ndarray: the passive force in Fourier space
        '''
        Pi_f = self.calc_passive_stress_f(Q)
        F_pf = []
        for j in range(self.dim):
            F_pf.append(np.sum(1j * self.k[i] *Pi_f[j][i] for i in range(self.dim)))
        return numpy.array(F_pf)

    def calc_passive_stress_f(self,Q):
        """
        Calculates the passive stress in fourier space

        Input:
            Q (numpy.narray) the order parameter that we use to find the stress.
        
        Output: 
            (numpy.narray) the passive stress in fourier space
        """
        if self.dim == 2:
            H = self.calc_molecular_field(Q)
            Antisym_QH = np.sum(self.get_sym_tl(Q,0,k)*self.get_sym_tl(H,k,1) -self.get_sym_tl(H,0,k)*self.get_sym_tl(Q,k,1) for k in range(self.dim))
            Ericksen = np.zeros((self.dim,self.xRes,self.yRes),dtype=np.complex128)
            Ericksen[0] = - self.K*np.sum(sp.fft.ifftn(1j*self.k[0]*sp.fft.fftn(self.get_sym_tl(Q,m,l)))*
                                              sp.fft.ifftn(1j * self.k[0] * sp.fft.fftn(self.get_sym_tl(Q,m,l)))
                                              for m in range(self.dim) for l in range(self.dim))
            Ericksen[1] = - self.K*np.sum(sp.fft.ifftn(1j*self.k[0]*sp.fft.fftn(self.get_sym_tl(Q,m,l)))*
                                              sp.fft.ifftn(1j * self.k[1] * sp.fft.fftn(self.get_sym_tl(Q,m,l)))
                                              for m in range(self.dim) for l in range(self.dim))

            stress = np.zeros((self.dim,self.dim,self.xRes,self.yRes),dtype=np.complex128)
            for i in range(self.dim):
                for j in range(self.dim):
                    stress[i][j] = self.get_sym_tl(Ericksen,i,j) + self.get_anti_sym(Antisym_QH,i,j)
            return sp.fft.fftn(stress, axes=(range(-self.dim, 0)) )

        elif self.dim == 3:
            # TODO optimise
            H = self.calc_molecular_field(Q)

            Antisym_QH = np.zeros((3, self.xRes, self.yRes,self.zRes), dtype=np.complex128)

            Antisym_QH[0] = np.sum(self.get_sym_tl(Q,0,k)*self.get_sym_tl(H,k,1) -self.get_sym_tl(H,0,k)*self.get_sym_tl(Q,k,1) for k in range(self.dim))
            Antisym_QH[1] = np.sum(
                self.get_sym_tl(Q, 0, k) * self.get_sym_tl(H, k, 2) - self.get_sym_tl(H, 0, k) * self.get_sym_tl(Q, k, 2) for k in
                range(self.dim))
            Antisym_QH[2] = np.sum(
                self.get_sym_tl(Q, 1, k) * self.get_sym_tl(H, k, 2) - self.get_sym_tl(H, 1, k) * self.get_sym_tl(Q, k, 2) for k in
                range(self.dim))

            Ericksen = np.zeros((5, self.xRes, self.yRes,self.zRes), dtype=np.complex128)
            Ericksen[0] = - self.K * np.sum(sp.fft.ifftn(1j * self.k[0] * sp.fft.fftn(self.get_sym_tl(Q, m, l))) *
                                            sp.fft.ifftn(1j * self.k[0] * sp.fft.fftn(self.get_sym_tl(Q, m, l)))
                                            for m in range(self.dim) for l in range(self.dim))
            Ericksen[1] = - self.K * np.sum(sp.fft.ifftn(1j * self.k[0] * sp.fft.fftn(self.get_sym_tl(Q, m, l))) *
                                            sp.fft.ifftn(1j * self.k[1] * sp.fft.fftn(self.get_sym_tl(Q, m, l)))
                                            for m in range(self.dim) for l in range(self.dim))
            Ericksen[2] = - self.K * np.sum(sp.fft.ifftn(1j * self.k[0] * sp.fft.fftn(self.get_sym_tl(Q, m, l))) *
                                            sp.fft.ifftn(1j * self.k[2] * sp.fft.fftn(self.get_sym_tl(Q, m, l)))
                                            for m in range(self.dim) for l in range(self.dim))
            Ericksen[3] = - self.K * np.sum(sp.fft.ifftn(1j * self.k[1] * sp.fft.fftn(self.get_sym_tl(Q, m, l))) *
                                            sp.fft.ifftn(1j * self.k[1] * sp.fft.fftn(self.get_sym_tl(Q, m, l)))
                                            for m in range(self.dim) for l in range(self.dim))
            Ericksen[4] = - self.K * np.sum(sp.fft.ifftn(1j * self.k[1] * sp.fft.fftn(self.get_sym_tl(Q, m, l))) *
                                            sp.fft.ifftn(1j * self.k[2] * sp.fft.fftn(self.get_sym_tl(Q, m, l)))
                                            for m in range(self.dim) for l in range(self.dim))

            stress = np.zeros((self.dim,self.dim,self.xRes,self.yRes,self.zRes))
            for i in range(self.dim):
                for j in range(self.dim):
                    stress[i][j] = self.get_sym_tl(Ericksen,i,j) + self.get_anti_sym(Antisym_QH,i,j)
            return sp.fft.fftn(stress, axes=(range(-self.dim, 0)) )

    def calc_trace_Q2(self,Q):

        if self.dim == 2:
            return 2*( Q[0]**2 + Q[1]**2)
        elif self.dim == 3:
            return 2 *(Q[0]**2 + Q[1]**2+ Q[2]**2+Q[3]**2 + Q[4]**2 + Q[0]*Q[3])


    def calc_molecular_field(self,Q):
        """
        Finds the molecular field (NB! need to be rewriten when C != 0)
        
        Input:
            Q (numpy.ndarray): The nematic tensor
        
        Output:
             (numpy.ndarray): The molecular field
        """

        Q2 =  self.calc_trace_Q2(Q)
        temp = -self.K * sp.fft.ifftn( self.k2* sp.fft.fftn(Q,axes=(range(-self.dim,0))),axes=(range(-self.dim,0)) )

        if self.dim == 2 or self.C == 0:
            return temp +self.A*self.B*Q -2*self.A*Q2*Q

        elif self.dim == 3:
            C_term = np.zeros((5,self.xRes,self.yRes,self.zRes))

            C_term[0] = self.C * (Q[0]**2 + Q[1]**2 +Q[2]**2 - 1/3 *Q2)
            C_term[1] = self.C * (Q[0]*Q[1] + Q[1]*Q[3] +Q[2]*Q[4])
            C_term[2] = self.C * (Q[1]*Q[4] -Q[2] * Q[3] )
            C_term[3] = self.C * (Q[1]**2 + Q[3]**2 +Q[4]**2 - 1/3 *Q2)
            C_term[4] = self.C * (Q[1]*Q[2]  -Q[4]*Q[0] )

            return temp +self.A*self.B*Q -2*self.A*Q2*Q +C_term


    def calc_pressure_f(self,F_af,F_pf):
        '''
        Calculates the pressure in Fourier space. The zero mode is set to zero
        
        Input:
            F_af (numpy.narray) the active force in Fourier space
            F_pf (numpy.narray) the passive force in Fourier space
        
        Output:
            (numpy.ndarray) the pressure
        '''
        p_af = np.sum(1j*self.k[i]*F_af[i] for i in range(self.dim))
        p_pf = np.sum(1j*self.k[i]*F_pf[i] for i in range(self.dim))
        return -(p_af + p_pf)/self.k2_press

    #TODO: This function needs a more descriptive name (Vidar 21.01.24)
    def calc_grad_p_f(self,p_f):
        """
        Caclulates the gradient of the pressure
        
        Input: 
            p_f (numpy.narray) the pressure in Fourier space
        
        Output:
             (numpy.ndarray) gradient of th epressure
        """
        grad_pf = []
        for i in range(self.dim):
            grad_pf.append(1j*self.k[i]*p_f)
        return np.array(grad_pf)

    def calc_vorticity_tensor(self):
        """
        Calculates the vorticity tensor

        Input:
            None

        Output:
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
        # TODO change this
        """
        Calculates the strainrate tensor

        Input:
            None

        Output:
            (numpy.ndarray) The strainrate
        """
        E_f = np.zeros_like(self.Q_f)
        for i in range(self.dim):
            for j in range(self.dim):
                E_f[i][j]= (1j*self.k[i]*self.u_f[j] +1j*self.k[j]*self.u_f[i])/2
        return E_f

#### Calculation of non-linear evolution terms
    def calc_nonlinear_evolution_function_f(self,Q,t):
        # TODO test and make sure that the passive stress works as intended (Jonas: 2023/11/14)
        """
        Calculates the non-linear evolution function for the nematic
        
        Input:
            Q (numpy.narray) the nematic orderparameter
        
        Output:
            (numpy.narray) the non-linear evolution function evaluated in Fourier space
        """
        if self.dim == 2:
            self.conf_u(Q)
            Q_f = sp.fft.fftn(Q,axes=range(-self.dim,0))
            N_f = self.calc_nonlinear_evolution_term_no_flow_f(Q,t)
            Omega =self.calc_vorticity_tensor()
            Antisym_Omega_Q = np.zeros_like(Q_f)

            Antisym_Omega_Q[0] = -2 *Q[1]*Omega
            Antisym_Omega_Q[1] = 2*Q[0]*Omega

            advectiv_deriv = - np.sum(self.u[k]* sp.fft.ifftn(1j*self.k[k] * Q_f,axes=(range(-self.dim,0)))for k in range(self.dim) )

            return sp.fft.fftn(Antisym_Omega_Q +advectiv_deriv, axes=range(-self.dim,0)) +N_f

        elif self.dim == 3:
            # TODO: check that antisym_omega is correct (18/01/24)
            self.conf_u(Q)
            Q_f = sp.fft.fftn(Q, axes=range(-self.dim, 0))
            N_f = self.calc_nonlinear_evolution_term_no_flow_f(Q,t)
            Omega = self.calc_vorticity_tensor()

            advectiv_deriv = - np.sum(
                self.u[k] * sp.fft.ifftn(1j * self.k[k] * Q_f, axes=(range(-self.dim, 0))) for k in range(self.dim))

            Antisym_Omega_Q = np.zeros_like(Q_f)

            Antisym_Omega_Q[0] = -2 * (Q[1]*Omega[0] + Q[2]*Omega[1] )
            Antisym_Omega_Q[1] = Omega[0]*(Q[0]-Q[3]) - Omega[2]*Q[2] -Omega[1]*Q[4]
            Antisym_Omega_Q[2] = Omega[1]*(2*Q[0] + Q[3] ) + Omega[1]*Q[1] - Omega[0]*Q[4]
            Antisym_Omega_Q[3] = 2*(Q[1]*Omega[0] -Q[4]*Omega[2])
            Antisym_Omega_Q[4] = Omega[2]*(2*Q[3]+Q[0]) + Omega[1]*Q[1] + Omega[0]*Q[2]

            return sp.fft.fftn(Antisym_Omega_Q + advectiv_deriv, axes=range(-self.dim, 0)) + N_f

        else:
            raise Exception("This dimension is not implemented at the moment")

    def calc_nonlinear_evolution_term_no_flow_f(self,Q,t):

        """
        Calculates the non-linear evolution function for the nematic without the flow field
        
        Input:
            Q (numpy.narray) the nematic orderparameter
        
        Output:
            (numpy.narray) the non-linear evolution function evaluated in Fourier space
        """
        Q2 = self.calc_trace_Q2(Q)
        if self.dim ==2 or self.C == 0:
            return -2*self.A*sp.fft.fftn(Q2 *Q,axes =(range(-self.dim,0)))/self.gamma

        elif self.dim == 3:
            C_term = np.zeros((5, self.xRes, self.yRes, self.zRes))

            C_term[0] = self.C * (Q[0] ** 2 + Q[1] ** 2 + Q[2] ** 2 - 1 / 3 * Q2)
            C_term[1] = self.C * (Q[0] * Q[1] + Q[1] * Q[3] + Q[2] * Q[4])
            C_term[2] = self.C * (Q[1] * Q[4] - Q[2] * Q[3])
            C_term[3] = self.C * (Q[1] ** 2 + Q[3] ** 2 + Q[4] ** 2 - 1 / 3 * Q2)
            C_term[4] = self.C * (Q[1] * Q[2] - Q[4] * Q[0])
            return -2*self.A*sp.fft.fftn(Q2 *Q,axes =(range(-self.dim,0)))/self.gamma + sp.fft.fftn(C_term ,axes =(range(-self.dim,0)))/self.gamma


##### evolvers
    def evolve_nematic(self, number_of_steps, method= 'ETD2RK'):
        '''
         Evolver for the nematic system
        
        Input:
            number_of_steps (int) the number of time steps that we are evolving the equation
            method (string, optional) the integration method we want to use. ETD2RK is sett as default
        
        Output:
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
        self.conf_u(self.Q)

    def evolve_nematic_no_flow(self,number_of_steps,method = 'ETD2RK'):
        ''' Evolver for the nematic system without the flow field
        
        Input:
            number_of_steps (int) the number of time steps that we are evolving the equation
            method (string, optional) the integration method we want to use. ETD2RK is sett as default
        
        Output:
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
    def calc_disclination_density_nematic(self):
        #TODO: See if this can be optimised
        """
        Calculates the defect density for the nematic. Note that in three dimension the defect density is a tensor

        Input:
            None

        Output:
            (numpy.narray) The defect density
        """
        if self.dim == 2:
            psi0 = np.sqrt(self.B)/2
            psi =[self.Q[0],self.Q[1]]
            return self.calc_defect_density(psi, psi0)

        elif self.dim == 3:
            D = np.zeros((self.dim,self.dim,self.xRes,self.yRes,self.zRes))
            term_trace = np.sum(np.real(sp.fft.ifftn(1j*self.k[k]* self.get_sym_tl(self.Q_f,k,a)))*
                                    np.real(sp.fft.ifftn(1j*self.k[l] * self.get_sym_tl(self.Q_f,l,a)))
                                    - np.real(sp.fft.ifftn(1j*self.k[k]* self.get_sym_tl(self.Q_f,l,a)))*
                                    np.real(sp.fft.ifftn(1j*self.k[l] * self.get_sym_tl(self.Q_f,k,a)))
                                    for k in range(self.dim) for a in range(self.dim) for l in range(self.dim))

            for gam in range(self.dim):
                for i in range(self.dim):
                    D[gam, i] = 2*np.sum(np.real(sp.fft.ifftn(1j*self.k[gam]* self.get_sym_tl(self.Q_f,k,a)))*
                                    np.real(sp.fft.ifftn(1j*self.k[k] * self.get_sym_tl(self.Q_f,i,a)))
                                    - np.real(sp.fft.ifftn(1j*self.k[gam]* self.get_sym_tl(self.Q_f,i,a)))*
                                    np.real(sp.fft.ifftn(1j*self.k[k] * self.get_sym_tl(self.Q_f,k,a)))
                                    for k in range(self.dim) for a in range(self.dim))
                    if gam == i:
                        D[gam,i] += term_trace
            return D
        else:
            raise Exception("Only two and three dimensions are currently supported")


    def calc_disclination_density_decoupled(self):
        if self.dim == 3:
            D = self.calc_disclination_density_nematic()
            S0 = self.calc_equilibrium_S()
            rho = D/(S0**2 *np.pi)
            omega = np.sqrt(np.sum(rho[i,j]*rho[i,j] for i in range(self.dim) for j in range(self.dim)) )

            DDT = np.zeros((self.xRes,self.yRes,self.zRes,self.dim,self.dim))
            DTD = np.zeros((self.xRes,self.yRes,self.zRes,self.dim,self.dim))

            for i in range(self.dim):
                for j in range(self.dim):
                    DDT[:,:,:,i,j] = np.sum(rho[i,k]*rho[j,k] for k in range(self.dim))
                    DTD[:, :, :, i, j] = np.sum(rho[k,i] * rho[ k,j] for k in range(self.dim))

            vals_1,vecs_1 =  numpy.linalg.eigh(DDT)
            vals_2, vecs_2 = numpy.linalg.eigh(DTD)

            Omega = np.transpose(vecs_1[:,:,:,:,2], (3,0,1,2))
            T = np.transpose(vecs_2[:,:,:,:,2], (3,0,1,2))
            trD = np.sum(D[i,i] for i in range(self.dim))
            return omega, Omega, T, trD


    def calc_dt_psi(self,Q_prev,delta_t):
        dt_Q = (self.Q -Q_prev)/delta_t
        return dt_Q[0] + 1j*dt_Q[1]

    def calc_equilibrium_S(self):
        '''
        Calculates the strength of nematic order S

        Input:
            None

        Output:
            (numpy.narray) equilibriums value of D
        '''
        if self.dim == 2:
            return  np.sqrt(self.B)

        elif self.dim == 3:
            S0 = 1 / 8 * self.C / self.A + 1 / 2 * np.sqrt(self.C ** 2 / (16 * self.A ** 2) + 3 * self.B)
            return S0

    def calc_order_and_director(self):
        """
        Finds the amount of order (S) and the director field (n)

        Input:
            None

        Output:
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
                    Q_eig[:,:,:,i,j] = self.get_sym_tl(self.Q,i,j)

            eigvals, eigvectors = numpy.linalg.eigh(Q_eig)
            S = 3/2 *eigvals[:,:,:,2]
            n = np.transpose(eigvectors[:,:,:,:,2], (3,0,1,2))

            return S, n

    def calc_vortex_velocity_field(self, dt_Q, psi=None):
        # TODO make 3D
        """
        Calculates the velocity field of the defect in two dimensions
        
        Input:
            dt_Q (numpy.narray) the time derivative of the order parameter
            psi (numpy.narray, optional) the order parameter on complex form
        
        Output:
            (numpy.narray) the velocity field
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

    def calc_disclination_nodes_nem(self, dt_Q=None,polarization = None,charge_tolerance=None):
        """
        Calculate the positions and charges of vortex nodes based on the defect density.
        
        Input:
            dt_Q (numpy.narray, optional): The time derivative of the order parameter. If not provided, the velocity of the vortex nodes will not be calculated.
        
        Output:
            list: A list of dictionaries representing the vortex nodes. Each dictionary contains the following keys:
                  - 'position_index': The position index of the vortex node in the defect density array.
                  - 'charge': The charge of the vortex node.
                  - 'position': The position of the vortex node as a list [x, y].
                  - 'velocity': The velocity of the vortex node as a list [vx, vy].
        """

        # Calculate defect density
        if self.dim == 2:
            psi = self.Q[0] + 1j*self.Q[1]
            rho = self.calc_disclination_density_nematic()

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

        elif self.dim == 3:
            #TODO Make sure that tangent vector is continous across the border
            omega, Omega, T, trD = self.calc_disclination_density_decoupled()

            vortex_nodes = self.calc_defect_nodes(omega,charge_tolerance=None)
            position_list = []

            for vortex in vortex_nodes:

                tangent_vector = np.array([T[i][vortex['position_index']] for i in range(3)])
                rotation_vector = np.array([Omega[i][vortex['position_index']] for i in range(3)])

                for i in range(len(position_list)):
                    pos = position_list[i]
                    if np.sqrt(sum( (vortex['position'][i] -pos[i])**2 for i in range(self.dim) )) <  5*self.a0:
                        tan_neight = vortex_nodes[i]['Tangent_vector']
                        if np.sum((tan_neight[j] -tangent_vector[j] )**2 -(tan_neight[j] + tangent_vector[j])**2 for j in range(self.dim)) > 0:

                            tangent_vector = -1* tangent_vector
                        break

                if np.sign(np.dot(tangent_vector, rotation_vector)) != np.sign(trD[vortex['position_index']]):
                    rotation_vector = -1*rotation_vector

                vortex['Tangent_vector'] = tangent_vector
                vortex['Rotation_vector'] = rotation_vector
                position_list.append(vortex['position'])



        return vortex_nodes

    def plot_disclination_nodes(self, vortex_nodes, ax=None):
        """
        Plots the discliation nodes in the given axes.

        Input:
            vortex_nodes (list): A list of dictionaries representing the vortex nodes. Each dictionary contains the following keys:
                                 - 'position': The position of the vortex node as a list [x, y].
                                 - 'charge': The charge of the vortex node.
                                 - 'velocity': The velocity of the vortex node as a list [vx, vy].
            ax (Axes, optional): The axes to plot the vortex nodes on. If not provided, a new subplot will be created.

        Output:
            matplotlib.axes.Axes: The axes on which the disclination nodes are plotted.
        """

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
            return ax
        elif self.dim == 3:
            # Plotting options
            quiver_scale = 2  # The scale of the quiver arrows

            if ax == None:
                ax = plt.gcf().add_subplot(111, projection='3d')
            x_coords = []
            y_coords = []
            z_coords = []

            tx = []
            ty = []
            tz = []

            Ox = []
            Oy = []
            Oz = []

            for vortex in vortex_nodes:
                x_coords.append(vortex['position'][0])
                y_coords.append(vortex['position'][1])
                z_coords.append(vortex['position'][2])

                tx.append(vortex['Tangent_vector'][0])
                ty.append(vortex['Tangent_vector'][1])
                tz.append(vortex['Tangent_vector'][2])

                Ox.append(vortex['Rotation_vector'][0])
                Oy.append(vortex['Rotation_vector'][1])
                Oz.append(vortex['Rotation_vector'][2])

            tx = np.array(tx)
            ty = np.array(ty)
            tz = np.array(tz)

            Ox = np.array(Ox)
            Oy = np.array(Oy)
            Oz = np.array(Oz)


            # ax.scatter(x_coords, y_coords, z_coords, marker='o', color='black')
            ax.quiver(x_coords, y_coords, z_coords, quiver_scale * tx, quiver_scale * ty, quiver_scale * tz,
                      color='blue')
            ax.quiver(x_coords, y_coords, z_coords, quiver_scale * Ox*0.75 , quiver_scale * Oy*0.75 ,
                      quiver_scale * Oz*0.75, color='green')

            ax.set_xlabel('$x/a_0$')
            ax.set_ylabel('$y/a_0$')
            ax.set_zlabel('$z/a_0$')

            ax.set_xlim([0, self.xmax - self.dx])
            ax.set_ylim([0, self.ymax - self.dy])
            ax.set_zlim([0, self.zmax - self.dz])

            ax.set_aspect('equal')
            ax.grid(True)

            return ax

    def plot_field_velocity_and_director(self, field, velocity, director, **kwargs):# ax=None, colorbar=True, colormap='viridis',
                                       #  cmax=None, cmin=None,
                                        # number_of_layers=1, hold=False):
        """
        Plot the field, velocity, and director in the given axes.

        Input:
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

        Output:
            ax (Axes): The axes with the plotted field, velocity, and director.

        Raises:
            Exception: If the plotting function is not yet configured for dimensions other than 2.

        Note: streamplot assumes xy indexing and not ij. I think it is suficient
        just transpose the matrices before putting them in
        """
        if field.dtype == bool:
         field = field.astype(float)

        # Check if the vector field is complex
        if np.iscomplexobj(field):
            print(
                 "\033[91mWarning: the provided field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.\033[0m")
            print('Max imaginary part: ', np.max(np.imag(field)))
        field = np.real(field)

        # Check if an axis object is provided
        fig = kwargs.get('fig', plt.gcf())
        ax = kwargs.get('ax', None)

        # Kewyord arguments
        colorbar = kwargs.get('colorbar', True)

        # Extend the field if not a complete array is given
        field = self.plot_tool_extend_field(field)

        if self.dim == 2:

            # Keyword arguments particular to the 2D case
            kwargs['grid'] = kwargs.get('grid', False)

            if ax == None:
                fig.clf()
                ax = plt.gcf().add_subplot(111)

            colormap = kwargs.get('colormap', 'viridis')

            if colormap == 'bluewhitered':
                colormap = tool_colormap_bluewhitered()

            elif colormap == 'sunburst':
                colormap = tool_colormap_sunburst()
            else:
                colormap = plt.get_cmap(colormap)

            X, Y = np.meshgrid(self.x, self.y, indexing='ij')
            pcm = ax.pcolormesh(X / self.a0, Y / self.a0, field, shading='gouraud', cmap=colormap)

            xlim = [self.xmin, self.xmax - self.dx]
            ylim = [self.ymin, self.ymax - self.dy]

            limits_provided = False
            if 'xlim' in kwargs:
                xlim = kwargs['xlim']
                limits_provided = True
            else:
                if 'xmin' in kwargs:
                    xlim[0] = kwargs['xmin']
                    limits_provided = True

                if 'xmax' in kwargs:
                    xlim[1] = kwargs['xmax']
                    limits_provided = True

            if 'ylim' in kwargs:
                ylim = kwargs['ylim']
                limits_provided = True
            else:
                if 'ymin' in kwargs:
                    ylim[0] = kwargs['ymin']
                    limits_provided = True

                if 'ymax' in kwargs:
                    ylim[1] = kwargs['ymax']
                    limits_provided = True
            if limits_provided:
                region_to_plot = np.zeros(self.dims).astype(bool)
                region_to_plot[(xlim[0] <= X) * (X <= xlim[1]) * (ylim[0] <= Y) * (Y <= ylim[1])] = True
                vlim = [np.min(field[region_to_plot]), np.max(field[region_to_plot])]

            else:
                vlim = [np.min(field), np.max(field)]

            # Set the value limitses
            if 'vlim' in kwargs:
                vlim = kwargs['vlim']
            else:
                if 'vmin' in kwargs:
                    vlim[0] = kwargs['vmin']
                if 'vmax' in kwargs:
                    vlim[1] = kwargs['vmax']

            if vlim[1] - vlim[0] < 1e-10:
                vlim = [vlim[0] - 0.05, vlim[1] + 0.05]

            pcm.set_clim(vmin=vlim[0], vmax=vlim[1])

            if 'vlim_symmetric' in kwargs:
                vlim_symmetric = kwargs['vlim_symmetric']
                if vlim_symmetric:
                    cmax = abs(field).max()
                    cmin = -cmax
                    pcm.set_clim(vmin=cmin, vmax=cmax)

            colorbar = kwargs.get('colorbar', True)

            if colorbar:
                cbar = plt.colorbar(pcm, ax=ax)



            ax.streamplot(X.T, Y.T, (velocity[0]).T, (velocity[1]).T, color='w')
            ax.quiver(X, Y, director[0], director[1], headwidth=0, scale=50)
            ax.quiver(X, Y, -director[0], -director[1], headwidth=0, scale=50)
            ax.set_aspect('equal')

        else:
            raise Exception("This plotting function not yet configured for other dimension")

        kwargs['ax'] = ax
        self.plot_tool_set_axis_properties(**kwargs)
        return fig, ax
