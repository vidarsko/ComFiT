import numpy
import numpy as np
from comfit.core.base_system import BaseSystem
import scipy as sp
from tqdm import tqdm
import matplotlib.pyplot as plt
from comfit.tool.tool_colormaps import tool_colormap_angle, tool_colormap_bluewhitered, tool_colormap_sunburst
from comfit.tool.tool_create_orthonormal_triad import tool_create_orthonormal_triad

from comfit.nematic_liquid_crystal.plot_disclination_nodes_matplotlib import plot_disclination_nodes_matplotlib
from comfit.nematic_liquid_crystal.plot_field_velocity_and_director_matplotlib import plot_field_velocity_and_director_matplotlib

from comfit.nematic_liquid_crystal.plot_disclination_nodes_plotly import plot_disclination_nodes_plotly
from comfit.nematic_liquid_crystal.plot_field_velocity_and_director_plotly import plot_field_velocity_and_director_plotly

from comfit.tool.tool_math_functions import levi_civita_symbol
from mpl_toolkits.mplot3d import axes3d
import matplotlib.cm as cm
from matplotlib.colors import Normalize

class NematicLiquidCrystal(BaseSystem):

    def __init__(self, dim, **kwargs):
        """Initializes a system to simulate a (active) nematic liquid crystal

        Args:
            dimension : The dimension of the system.
            kwargs : dict, optional
            Optional keyword arguments to set additional parameters. See
            https://comfitlib.com/ClassNematicLiquidCrystal/

        Returns:
        - NematicLiquidCrystal object
            The system object representing the nematic simulation.

        Example:
            nematic = NematicLiquidCrystal(2, xRes=100, yRes = 100, alpha=-0.5)
            Creates a nematic liquid crystal with 2 dimensions and a spatial resolution of 100.
            The activity alpha is set to -0.5.
        """
        super().__init__(dim, **kwargs)

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
        """Output a string representation of the system.

        Input:
            None

        Returns:
            A string representation of the system.
        """
        return "NematicLiquidCrystal"

    # Initial condition
    def conf_initial_condition_ordered(self, noise_strength=0.01):
        """Configures the system with the nematogens pointing in the x-direction in 2D and in the z-direction in 3D
        with some random noise in the angle.
        
        Args:
             noise_strength: A meshure for how much noise to put in the angle (float)
        
        Returns:
            Initialises self.Q and self.Q_f
        
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
        """Sets the initial condition for a disclination dipole configuration in a 2-dimensional system.

        Args:
            None

        Returns:
            Configures self.Q and self.Q_f with a disclination dipole configuration.

        Raises:
            Exception: If the dimension of the system is not 2.
        """
        if not (self.dim == 2):
            raise Exception("The dimension of the system must be 2 for a disclination dipole configuration.")

        if self.Q is None:
            self.conf_initial_condition_ordered(noise_strength=0)

        if dipole_vector is None:
            dipole_vector = [(self.xmax -self.xmin)/ 3, 0]

        if dipole_position is None:
            dipole_position = self.rmid

        psi = (self.Q[0] + 1j*self.Q[1]) * np.exp(1j * self.calc_angle_field_vortex_dipole(dipole_vector, dipole_position))
        self.Q[0] = np.real(psi)
        self.Q[1] = np.imag(psi)
        self.Q_f = sp.fft.fft2(self.Q)

    def conf_initial_disclination_lines(self, position1=None,position2 = None):
        """Sets the initial condition for a disclination line in a 3-dimensional system.

        The dislocation is parralell to the z-axis

        Args:
            position1 (list): the position of the first dislocation. Only the position in the xy plane is used

        Returns:
            Sets the value of self.Q and self.Q_f
        """
        if not (self.dim == 3):
            raise Exception("The dimension of the system must be 3 for a disclination line configuration.")

        if position1 is None:
            position1 = [self.xmid+(self.xmax-self.xmin)/3, self.ymid]
        if position2 is None:
            position2 = [self.xmid - (self.xmax-self.xmin) / 3, self.ymid]


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




    def conf_active_channel(self,width = None,d=7):
        """Configures the activity to zero everywhere exept for inside a channel of width "width"
        
        Args:
            width: width of the channel (float)
            d: width of interface (float, optional)
        
        Returns:
            Updates the activity to the channel configuration.
        """
        if self.dim ==2:
            X, Y = np.meshgrid(self.x, self.y, indexing='ij')
            alpha_0 = self.alpha
            if width is None:
                width = (self.xmax-self.xmin)/3

            self.alpha = alpha_0*(1- 1 / 2 * (2 + np.tanh((X - self.xmid - width/2) / d) - np.tanh((X - self.xmid + width/2) / d)))
        else:
            raise Exception("The active channel is only permitted in two dimensions")

    ## calculations related to the flow field
    def conf_velocity(self,Q):
        '''
        Updates the velocity and its fourier transform given a nematic field Q.
        
        Args:
            (numpy.narray) the Q tensor
        
        Returns:
            (numpy.narray) velocity
        '''
        F_af = self.calc_active_force_f(Q)
        F_pf = self.calc_passive_force_f(Q)
        p_f = self.calc_pressure_f(F_af,F_pf)
        grad_pf = self.calc_gradient_pressure_f(p_f)
        self.u_f = (F_af + F_pf-grad_pf )/ (self.Gamma +self.eta*self.k2_press)

        self.u = np.real(sp.fft.ifftn(self.u_f, axes=(range(-self.dim, 0))))

    def calc_active_force_f(self,Q):
        '''Function that calculates the activ force in Fourier space.
        
        Args:
            Q: the order parameter that we use to find the force.  (numpy.narray) 
        
        Returns:
            The active force in Fourier space (numpy.narray) 
        '''
        F_af = []
        for j in range(self.dim):
            F_af.append(np.sum(1j*self.k[i]*sp.fft.fftn(self.alpha *self.get_sym_tl(Q,j,i),axes=(range(-self.dim, 0)) ) for i in range(self.dim)))
        return np.array(F_af)

    def calc_passive_force_f(self,Q):
        '''Calculates the passive force in Fourier space
        
        Args:
            Q: the order parameter that we use to find the force. (numpy.narray)
        
        Returns: 
            The passive force in Fourier space numpy.ndarray: 
        '''
        Pi_f = self.calc_passive_stress_f(Q)
        F_pf = []
        for j in range(self.dim):
            F_pf.append(np.sum(1j * self.k[i] *Pi_f[j][i] for i in range(self.dim)))
        return numpy.array(F_pf)

    def calc_passive_stress_f(self,Q):
        """Calculates the passive stress in fourier space

        Args:
            Q: the order parameter that we use to find the stress.
        
        Returns: 
            The passive stress in fourier space (numpy.narray) 
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
        """Finds the molecular field 
        
        Args:
            Q (numpy.ndarray): The nematic tensor
        
        Returns:
            The molecular field (numpy.ndarray)
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
        '''Calculates the pressure in Fourier space. The zero mode is set to zero
        
        Args:
            F_af: the active force in Fourier space  (numpy.narray)
            F_pf:the passive force in Fourier space  (numpy.narray) 
        
        Returns: 
            The pressure (numpy.ndarray) 
        '''
        p_af = np.sum(1j*self.k[i]*F_af[i] for i in range(self.dim))
        p_pf = np.sum(1j*self.k[i]*F_pf[i] for i in range(self.dim))
        return -(p_af + p_pf)/self.k2_press

    
    def calc_gradient_pressure_f(self,p_f):
        """Caclulates the gradient of the pressure
        
        Args: 
            p_f: the pressure in Fourier space  (numpy.narray)
        
        Returns:
            Gradient of the pressure (numpy.ndarray) 
        """
        grad_pf = []
        for i in range(self.dim):
            grad_pf.append(1j*self.k[i]*p_f)
        return np.array(grad_pf)

    def calc_vorticity_tensor(self):
        """Calculates the vorticity tensor

        Args:
            None

        Returns:
            The vorticity tensor (numpy.ndarray) 
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
        """
        Calculates the strainrate tensor

        Args:
            None

        Returns: 
            The strainrate (numpy.ndarray) 
        """
        trace_u = np.sum(1j*self.k[i]*self.u_f[i] for i in range(self.dim))
        if self.dim == 2:
            E_f = np.zeros((2,self.xRes,self.yRes),dtype=np.complex128)
            E_f[0] = (1j*self.k[0] *self.u_f[0]) - trace_u/2
            E_f[1] = (1j*self.k[0] *self.u_f[1] +1j*self.k[1] *self.u_f[0])/2
        elif self.dim == 3:
            E_f = np.zeros((5, self.xRes, self.yRes,self.zRes), dtype=np.complex128)
            E_f[0] = (1j*self.k[0] *self.u_f[0]) - trace_u/3
            E_f[1] = (1j*self.k[0] *self.u_f[1] +1j*self.k[1] *self.u_f[0])/2
            E_f[2] = (1j * self.k[0] * self.u_f[2] + 1j * self.k[2] * self.u_f[0]) / 2
            E_f[3] =  (1j*self.k[1] *self.u_f[1]) - trace_u/3
            E_f[4] = (1j * self.k[1] * self.u_f[2] + 1j * self.k[2] * self.u_f[1]) / 2
        return E_f

    ## Calculation of non-linear evolution terms
    def calc_nonlinear_evolution_function_f(self,Q,t):
        """Calculates the non-linear evolution function for the nematic
        
        Args:
            Q: the nematic orderparameter (numpy.narray) 
        
        Returns:
            
            The non-linear evolution function evaluated in Fourier space (numpy.narray) 
        """
        if self.dim == 2:
            self.conf_velocity(Q)
            Q_f = sp.fft.fftn(Q,axes=range(-self.dim,0))
            N_f = self.calc_nonlinear_evolution_term_no_flow_f(Q,t)
            Omega =self.calc_vorticity_tensor()
            Antisym_Omega_Q = np.zeros_like(Q_f)

            Antisym_Omega_Q[0] = -2 *Q[1]*Omega
            Antisym_Omega_Q[1] = 2*Q[0]*Omega

            advectiv_deriv = - np.sum(self.u[k]* sp.fft.ifftn(1j*self.k[k] * Q_f,axes=(range(-self.dim,0)))for k in range(self.dim) )

            return sp.fft.fftn(Antisym_Omega_Q +advectiv_deriv, axes=range(-self.dim,0)) +N_f

        elif self.dim == 3:
            self.conf_velocity(Q)
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
        """Calculates the non-linear evolution function for the nematic without the flow field
        
        Args:
            Q: the nematic orderparameter  (numpy.narray)
        
        Returns:
            The non-linear evolution function evaluated in Fourier space (numpy.narray) 
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
        '''Evolves the nematic system
        
        Args:
            number_of_steps: the number of time steps that we are evolving the equation  (int)
            method: the integration method we want to use. ETD2RK is sett as default  (string)
        
        Returns:
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
        self.conf_velocity(self.Q)

    def evolve_nematic_no_flow(self,number_of_steps,method = 'ETD2RK'):
        ''' Evolves the nematic system without the flow field
        
        Args:
            number_of_steps: the number of time steps that we are evolving the equation  (int)
            method: the integration method we want to use. ETD2RK is sett as default  (string)
        
        Returns:
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


    ## Disclination tracking
    def calc_disclination_density_nematic(self):
        """
        Calculates the disclination density for the nematic. Note that in three dimension the disclination density is a tensor

        Args:
            None

        Returns:
            The disclination density (numpy.narray) 
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
            S0 = self.calc_equilibrium_S()
            return D/(S0**2 *np.pi)
        else:
            raise Exception("Only two and three dimensions are currently supported")


    def calc_disclination_density_decoupled(self):
        """Calculates the decoupled disclination density

        Args:
            None
        
        Returns:
            The disclination density (numpy.narray)
        """ 
        if self.dim == 3:
            rho = self.calc_disclination_density_nematic()


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
            trRho = np.sum(rho[i,i] for i in range(self.dim))
            return omega, Omega, T, trRho


    def calc_dt_psi(self,Q_prev,delta_t):
        """Calculates the time derivative of the order parameter as a complex field

        Args:
            Q_prev: the order parameter at the previous time step (numpy.narray)
            delta_t: the time step (float)
        
        Returns:
            The time derivative of the order parameter (numpy.narray)
        """

        dt_Q = (self.Q -Q_prev)/delta_t
        return dt_Q[0] + 1j*dt_Q[1]

    def calc_equilibrium_S(self):
        '''Calculates the strength of nematic order S

        Args:
            None

        Returns:
            equilibrium value of S (numpy.narray)
        '''
        if self.dim == 2:
            return  np.sqrt(self.B)

        elif self.dim == 3:
            S0 = 1 / 8 * self.C / self.A + 1 / 2 * np.sqrt(self.C ** 2 / (16 * self.A ** 2) + 3 * self.B)
            return S0

    def calc_order_and_director(self):
        """Calculates the amount of order (S) and the director field (n)

        Args:
            None

        Returns:
            Tuple consisting of 
                - Amount of order (scalar field) 
                - the director field (vector field)
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

    def calc_disclination_velocity_field(self, dt_Q):
        """
        Calculates the velocity field of the disclination in two dimensions
        
        Args:
            dt_Q: the time derivative of the order parameter (numpy.narray)

        Returns:
            The velocity field (numpy.narray) 
        """
        if self.dim ==2:

            psi = self.Q[0] +1j*self.Q[1]

            dt_psi = dt_Q[0] + 1j * dt_Q[1]

            return self.calc_defect_velocity_field([np.real(psi), np.imag(psi)],
                                                   [np.real(dt_psi), np.imag(dt_psi)])


    def calc_disclination_polarization_field(self):
        """Calculates the polarization field of the disclination in two dimensions

        Args:
            None
        
        Returns:
            The polarization field (numpy.narray)
        """
        ex = np.real(sp.fft.ifftn(1j*self.k[0]*self.Q_f[0] + 1j*self.k[1]*self.Q_f[1]))
        ey = np.real(sp.fft.ifftn(1j * self.k[0] * self.Q_f[1] - 1j * self.k[1] * self.Q_f[0]))
        return np.array([ex,ey])

    def calc_disclination_nodes_nem(self, dt_Q=None,polarization = None,charge_tolerance=None):
        """Calculates the positions and charges of disclination nodes based on the disclination density.
        
        Args:
            dt_Q: The time derivative of the order parameter. If not provided, the velocity of the disclination nodes will not be calculated. (numpy.narray, optional)
        
        Returns:
            A list of dictionaries representing the disclination nodes. Each dictionary contains the following keys:
                  - 'position_index': The position index of the disclination node in the disclination density array.
                  - 'charge': The charge of the disclination node.
                  - 'position': The position of the disclination node as a list [x, y].
                  - 'velocity': The velocity of the disclination node as a list [vx, vy].
        """

        # Calculate disclination density
        if self.dim == 2:

            rho = self.calc_disclination_density_nematic()

            if dt_Q is not None:
                velocity_field = self.calc_disclination_velocity_field(dt_Q)


            disclination_nodes = self.calc_defect_nodes(np.abs(rho))
            for disclination in disclination_nodes:
                disclination['charge'] = np.sign(rho[disclination['position_index']])

                if dt_Q is not None:
                    disclination['velocity'] = [velocity_field[0][disclination['position_index']],
                                          velocity_field[1][disclination['position_index']]]
                else:
                    disclination['velocity'] = [float('nan'), float('nan')]

                if disclination['charge'] >0 and polarization is not None:
                    disclination['polarization'] = [polarization[0][disclination['position_index']]
                        ,polarization[1][disclination['position_index']]]
                else:
                    disclination['polarization'] = [float('nan'), float('nan')]

        elif self.dim == 3:
            omega, Omega, T, trD = self.calc_disclination_density_decoupled()
            S0 = self.calc_equilibrium_S()
            disclination_nodes = self.calc_defect_nodes(omega,charge_tolerance=S0/np.pi)
            position_list = []

            for disclination in disclination_nodes:

                tangent_vector = np.array([T[i][disclination['position_index']] for i in range(3)])
                rotation_vector = np.array([Omega[i][disclination['position_index']] for i in range(3)])

                for i in range(len(position_list)):
                    pos = position_list[i]
                    if np.sqrt(sum( (disclination['position'][i] -pos[i])**2 for i in range(self.dim) )) <  5*self.a0:
                        tan_neight = disclination_nodes[i]['Tangent_vector']
                        if np.sum((tan_neight[j] -tangent_vector[j] )**2 -(tan_neight[j] + tangent_vector[j])**2 for j in range(self.dim)) > 0:

                            tangent_vector = -1* tangent_vector
                        break

                if np.sign(np.dot(tangent_vector, rotation_vector)) != np.sign(trD[disclination['position_index']]):
                    rotation_vector = -1*rotation_vector

                disclination['Tangent_vector'] = tangent_vector
                disclination['Rotation_vector'] = rotation_vector
                position_list.append(disclination['position'])


        return disclination_nodes

    def plot_disclination_nodes(self, disclination_nodes, **kwargs):
        if self.plot_lib == 'plotly':
            return plot_disclination_nodes_plotly(self, disclination_nodes, **kwargs)
        elif self.plot_lib == 'matplotlib':
            return plot_disclination_nodes_matplotlib(self, disclination_nodes, **kwargs)

    def plot_field_velocity_and_director(self, field, velocity, director, **kwargs):
        if self.plot_lib == 'plotly':
            return plot_field_velocity_and_director_plotly(self, field, velocity, director, **kwargs)
        elif self.plot_lib == 'matplotlib':
            return plot_field_velocity_and_director_matplotlib(self, field, velocity, director, **kwargs)
