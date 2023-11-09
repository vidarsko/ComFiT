import numpy as np
import matplotlib.pyplot as plt
from comfit.core.base_system import BaseSystem


class BEC(BaseSystem):
    def __init__(self, dimension, **kwargs):
        """
        Initializes a system to simulate a Bose-Einstein Condensate using the Gross-Pitaevskii equation.

        Parameters:
        - dimension : int
            The dimension of the system.
        - x_resolution : int
            The resolution along the x-axis.
        - kwargs : dict, optional
            Optional keyword arguments to set additional parameters.

        Returns:
        - BEC object
            The system object representing the BEC simulation.

        Example:
        bec = BEC(3, 100, dissipation_factor=0.5)
        Creates a BEC system with 3 dimensions and an x-resolution of 100. The dissipation factor is set to 0.5.
        """

        # First initialize the BaseSystem class
        super().__init__(dimension, **kwargs)

        # Type of the system
        self.psi = None
        self.psi_f = None
        self.type = 'BEC'

        # Default simulation parameters
        self.gamma = 0.1 if 'gamma' not in kwargs else kwargs['gamma']    # Dissipation (gamma)
        self.V_ext = 0  # External potential

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

    # Initial conditions
    def set_initial_condition_disordered(self, noise_strength=0.01):
        """
        Sets disordered initial condition for the BEC with some thermal flcutiations

        :param noise_strength:
        :return:
        """

        if self.dim==1:
            self.psi = np.random.rand(self.xRes) - 0.5 + \
                       1j * np.random.rand(self.xRes) - 0.5j

        elif self.dim==2:
            self.psi = np.random.rand(self.xRes, self.yRes)-0.5 + \
                       1j*np.random.rand(self.xRes, self.yRes)-0.5j

        elif self.dim==3:
            self.psi = np.random.rand(self.xRes, self.yRes,self.zRes) - 0.5 + \
                       1j * np.random.rand(self.xRes, self.yRes,self.zRes) - 0.5j
        else:
            raise Exception("Code for this dimension has not yet been implemented.")

        self.psi = noise_strength*self.psi
        self.psi_f = np.fft.fftn(self.psi)

    def set_harmonic_potential(self,R_tf):
        """
        Set the external potential to a harmonic trap with R_tf being the thomas fermi radius
        :param R_tf:
        :return:
        """
        trapping_strength = 1 / (R_tf ** 2)
        if self.dim ==1:
            return  trapping_strength*(self.x -self.xmid )**2
        if self.dim == 2:
            return trapping_strength*(((self.x-self.xmid)**2).reshape(self.xRes, 1)
                                         +((self.y-self.ymid)**2).reshape(1, self.yRes) )
        if self.dim == 3:
            return trapping_strength * (((self.x - self.xmid) ** 2).reshape(self.xRes, 1,1)
                                           + ((self.y - self.ymid) ** 2).reshape(1, self.yRes,1)
                                           +((self.z - self.zmid) ** 2).reshape(1, 1,self.zRes) )

    def set_initial_condition_Thomas_Fermi(self):
        """
        Finds the Thomas_Fermi ground state of a given potential
        must be precided by an energy relaxation to find the true ground state
        :return:
        """
        self.psi = np.emath.sqrt(1-self.V_ext)
        self.psi[self.V_ext > 1] = 0
        self.psi_f = np.fft.fftn(self.psi)

    def set_spatialy_varying_gamma(self,d=7, wx=0,wy=0,wz=0):
        '''
        this function sets a dissipative frame around the computational domain
        :param d: length of the interface between the low gamma and high gamma regions
        :param wx: distance fom center to the frame in x-direction
        :param wy:     -- " --                         y-direction
        :param wz:     -- " --                         z-direction
        :return:
        '''
        if self.dim == 2:
            X,Y =  np.meshgrid(self.x, self.y, indexing='ij')
            gammax = self.gamma +1/2 * (2 + np.tanh((X-self.xmid-wx)/d) - np.tanh((X-self.xmid+wx)/d) )
            gammay = self.gamma + 1 / 2 * (2 + np.tanh((Y-self.ymid - wy) / d) - np.tanh((Y-self.ymid + wy) / d))
            self.gamma = np.real(np.maximum(gammax,gammay))
        elif self.dim == 3:
            X, Y, Z = np.meshgrid(self.x, self.y,self.z, indexing='ij')
            gammax = self.gamma + 1 / 2 * (2 + np.tanh((X-self.xmid - wx) / d) - np.tanh((X-self.xmid + wx) / d))
            gammay = self.gamma + 1 / 2 * (2 + np.tanh((Y-self.ymid - wy) / d) - np.tanh((Y-self.ymid + wy) / d))
            gammaz = self.gamma + 1 / 2 * (2 + np.tanh((Z-self.zmid - wz) / d) - np.tanh((Z-self.zmid + wz) / d))
            self.gamma = np.real(np.maximum(gammax, gammay,gammaz))
        else:
            raise Exception("This feature is not yet available for the given dimension.")

    def gaussian_stirring_potential(self,size,strength,position):

        if self.dim ==1:
            return strength*np.exp( -(self.x -position[0])**2/size**2 )

        elif self.dim == 2:
            return strength*np.exp(-(((self.x-position[0])**2).reshape(self.xRes, 1)
                                         + ((self.y-position[1])**2).reshape(1, self.yRes)) /(size**2) )
        elif self.dim == 3:
            return strength* np.exp(-(((self.x - position[0]) ** 2).reshape(self.xRes, 1,1)
                                           + ((self.y - position[1]) ** 2).reshape(1, self.yRes,1)
                                           +((self.z - position[2]) ** 2).reshape(1, 1,self.zRes))/(size**2 ))

    #Time evolution
    def evolve_dGPE(self,number_of_steps):
        k2 = self.calc_k2()
        omega_f =   (1j + self.gamma) * (1 - 1 / 2 * k2)
        integrating_factors_f = self.calc_evolution_integrating_factors_ETDRK2(omega_f)

        for n in range(number_of_steps):
            self.psi, self.psi_f = self.evolve_ETDRK2_loop(integrating_factors_f,self.calc_nonlinear_evolution_term_f,
                                                           self.psi, self.psi_f)

    def evolve_dGPE_ETDRK4(self,number_of_steps):
        k2 = self.calc_k2()
        omega_f = (1j + self.gamma) * (1 - 1 / 2 * k2)
        integrating_factors_f = self.calc_evolution_integrating_factors_ETDRK4(omega_f)
        for n in range(number_of_steps):
            self.psi, self.psi_f = self.evolve_ETDRK4_loop(integrating_factors_f, self.calc_nonlinear_evolution_term_f,
                                                           self.psi, self.psi_f)




    def evolve_time_dependent(self,number_of_steps,V_t):
        """ evolver for time dependent potential, takes inn a function V_t
         that is the time dependent potential.
          Gives the user the freedom to use whatever potential they want"""
        if not(hasattr(self,'t')):
            self.t = 0
        k2 = self.calc_k2()

        omega_f = (1j + self.gamma) * (1 - 1 / 2 * k2)
        integrating_factors_f = self.calc_evolution_integrating_factors_ETDRK2(omega_f)
        for n in range(number_of_steps):
            self.psi, self.psi_f = self.evolve_ETDRK2_loop_timedep(integrating_factors_f,self.calc_nonlinear_evolution_term_timedep_f,
                                                                    V_t, self.psi, self.psi_f)

    def evolve_time_dependent_ETDRK4(self,number_of_steps,V_t):
        if not (hasattr(self, 't')):
            self.t = 0
        k2 = self.calc_k2()
        omega_f = (1j + self.gamma) * (1 - 1 / 2 * k2)
        integrating_factors_f = self.calc_evolution_integrating_factors_ETDRK4(omega_f)
        for n in range(number_of_steps):
            self.psi, self.psi_f = self.evolve_ETDRK4_loop_timedep(integrating_factors_f,
                                                                   self.calc_nonlinear_evolution_term_timedep_f,
                                                                   V_t, self.psi, self.psi_f)

    def evolve_relax_BEC(self,number_of_steps):
        """
        Evolves the BEC in 'imaginary time' to reach a stable (low free energy state).
        :param number_of_steps:
        :return:
        """
        gamma0 = self.gamma

        self.gamma=1-1j

        self.evolve_dGPE_ETDRK4(number_of_steps)

        self.gamma=gamma0

    def evolve_comoving_dGPE(self, number_of_steps, velx):
        #TODO This one is unstable in 3D (2D is fine) when using ETDRK4 (no problem with ETDRK2)
        # i am not sure why, and this needs to be figured out. The ETDRK4 do not have stability
        # problems in when stiring in the normal frame in 3D
        k2 = self.calc_k2()

        omega_f = (1j) * (1 - 1 / 2 * k2) + velx * self.dif[0]

        integrating_factors_f = self.calc_evolution_integrating_factors_ETDRK2(omega_f)

        for n in range(number_of_steps):
            self.psi, self.psi_f = self.evolve_ETDRK2_loop(integrating_factors_f, self.calc_nonlinear_evolution_term_comoving_f,
                                                           self.psi, self.psi_f)


            #Calculation functions


    def calc_nonlinear_evolution_term_f(self,psi):
        psi2 = np.abs(psi)**2
        return np.fft.fftn((1j+self.gamma)*(-self.V_ext-psi2)*psi)

    def calc_nonlinear_evolution_term_timedep_f(self,psi,V_t):
        self.V_ext = V_t(self.t)
        psi2 = np.abs(psi) ** 2
        return np.fft.fftn((1j + self.gamma) * (-self.V_ext - psi2) * psi)

    def calc_nonlinear_evolution_term_comoving_f(self,psi):
        psi2 = np.abs(psi)**2
        term1 = np.fft.fftn(-(1j+self.gamma)*(self.V_ext+psi2)*psi)
        term2 = np.fft.fftn(self.gamma*psi)
        term3 = 0.5*np.fft.fftn(self.gamma*np.fft.ifftn(-self.calc_k2()*np.fft.fftn(psi)))
        return (term1 + term2 + term3)

    def calc_vortex_density(self):
        return self.calc_defect_density([np.real(self.psi),np.imag(self.psi)])

    def calc_vortex_density_singular(self):
        return self.calc_defect_density([np.real(self.psi),np.imag(self.psi)])




    def calc_vortex_nodes(self):

        charge_tolerance = 0.2

        #Calculate defect density
        rho = self.calc_vortex_density()

        #Calculate the point where defect density is largest
        rho_max_index = np.unravel_index(np.argmax(np.abs(rho)),rho.shape)

        #Integrate the defect density around this point (i.e. in a ball around)
        charge,ball = self.calc_integrate_field(rho,index=rho_max_index,radius=1)

        #self.plot_field(rho)
        #plt.show()

        vortex_nodes = []

        X, Y = np.meshgrid(self.x, self.y, indexing='ij')


        while abs(charge)>charge_tolerance:
            vortex = {}
            vortex['position_index'] = rho_max_index
            vortex['charge'] = np.sign(charge)*np.ceil(np.abs(charge))
            x = np.sum(ball*abs(rho)*X)/np.sum(ball*abs(rho))
            y = np.sum(ball*abs(rho)*Y)/np.sum(ball*abs(rho))
            vortex['position'] = [x,y]

            vortex_nodes.append(vortex)

            rho[ball]=0
            rho_max_index = np.unravel_index(np.argmax(np.abs(rho)), rho.shape)

            charge, ball = self.calc_integrate_field(rho, index=rho_max_index, radius=1)



        return vortex_nodes


    #Plot functions
    def plot_vortices(self,vortex_nodes,ax=None):

        if ax == None:
            ax = plt.gcf().add_subplot(111)

        x_coords_pos = []
        y_coords_pos = []

        x_coords_neg = []
        y_coords_neg = []

        for vortex in vortex_nodes:

            if vortex['charge'] > 0:
                x_coords_pos.append(vortex['position'][0])
                y_coords_pos.append(vortex['position'][1])
            else:
                x_coords_neg.append(vortex['position'][0])
                y_coords_neg.append(vortex['position'][1])

        print(x_coords_pos,y_coords_pos)
        print(x_coords_neg,y_coords_neg)
        ax.scatter(x_coords_pos,y_coords_pos, marker='+',color='red')
        ax.scatter(x_coords_neg,y_coords_neg, marker='*',color='blue')
        ax.set_aspect('equal')

        

