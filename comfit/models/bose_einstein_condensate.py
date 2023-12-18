import numpy as np
import matplotlib.pyplot as plt
from comfit.core.base_system import BaseSystem
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D  # for 3D plotting


class BoseEinsteinCondensate(BaseSystem):
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
        - BoseEinsteinCondensate object
            The system object representing the BoseEinsteinCondensate simulation.

        Example:
        bec = BoseEinsteinCondensate(3, 100, dissipation_factor=0.5)
        Creates a BoseEinsteinCondensate system with 3 dimensions and an x-resolution of 100. The dissipation factor is set to 0.5.
        """

        # First initialize the BaseSystem class
        super().__init__(dimension, **kwargs)

        # Type of the system
        self.psi = None
        self.psi_f = None
        self.type = 'BoseEinsteinCondensate'

        # Default simulation parameters
        self.gamma = 0.1 if 'gamma' not in kwargs else kwargs['gamma']  # Dissipation (gamma)

        self.V0 = 0
        self.V_ext = lambda: self.V0  # External potential, note this is a function

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

    # SETTING FUNCTIONS
    def conf_initial_condition_disordered(self, noise_strength=0.01):
        """
        Sets disordered initial condition for the BoseEinsteinCondensate with some thermal flcutiations

        :param noise_strength:
        :return:
        """

        if self.dim == 1:
            self.psi = np.random.rand(self.xRes) - 0.5 + \
                       1j * np.random.rand(self.xRes) - 0.5j

        elif self.dim == 2:
            self.psi = np.random.rand(self.xRes, self.yRes) - 0.5 + \
                       1j * np.random.rand(self.xRes, self.yRes) - 0.5j

        elif self.dim == 3:
            self.psi = np.random.rand(self.xRes, self.yRes, self.zRes) - 0.5 + \
                       1j * np.random.rand(self.xRes, self.yRes, self.zRes) - 0.5j
        else:
            raise Exception("Code for this dimension has not yet been implemented.")

        self.psi = noise_strength * self.psi
        self.psi_f = np.fft.fftn(self.psi)


    def conf_time_dependent_potential(self, Func):
        """
        Set the potential to the function Func. Func has to use self.dt as the time variabel.
        Args:
            Func (function): the time-dependent potetnial that we want
        returns:
            Set V_ext to Func
        """
        self.V_ext = Func

    def conf_initial_condition_Thomas_Fermi(self):
        """
        Finds the Thomas_Fermi ground state.
        Must be precided by an energy relaxation to find the true ground state

        Returns: sets the value of self.psi and self.psi_f
        """
        V_0 = np.zeros(self.dims) + self.V_ext()
        self.psi = np.emath.sqrt(1 - V_0)

        self.psi[V_0 > 1] = 0
        self.psi_f = np.fft.fftn(self.psi)

    # CONFIGURATION FUNCTIONS
    def conf_insert_vortex(self, charge=1, position=None):
        """
        Sets the initial condition for a vortex dipole
        Returns:
        Modifies the value of self.psi and self.psi_f
        """
        if not (self.dim == 2):
            raise Exception("The dimension of the system must be 2 for a single point vortex.")

        if position is None:
            position = [self.xmid, self.ymid]

        if self.psi is None:
            self.conf_initial_condition_Thomas_Fermi()

            # TODO: maybe this needs to be formulated in terms of model parameters (Vidar 16.11.23)
            #  Answer: Homogeneous ground-state is now replaced by the Thomas-Fermi ground-state (Jonas 21.11.23 )

        self.psi = self.psi * np.exp(1j * self.calc_angle_field_single_vortex(position, charge=charge))
        self.psi_f = np.fft.fftn(self.psi)

    def conf_insert_vortex_dipole(self, dipole_vector=None, dipole_position=None):
        """
        Sets the initial condition for a vortex dipole configuration in a 2-dimensional system.

        Parameters:
            None

        Returns:
            None

        Raises:
            Exception: If the dimension of the system is not 2.
        """
        if not (self.dim == 2):
            raise Exception("The dimension of the system must be 2 for a vortex dipole configuration.")

        if self.psi is None:
            self.conf_initial_condition_Thomas_Fermi()

        if dipole_vector is None:
            dipole_vector = [self.xmax / 3, 0]

        if dipole_position is None:
            dipole_position = self.rmid

        self.psi = self.psi * np.exp(1j * self.calc_angle_field_vortex_dipole(dipole_vector, dipole_position))
        self.psi_f = np.fft.fftn(self.psi)

    def conf_insert_vortex_ring(self, position=None, radius=None, normal_vector=[0, 0, 1]):
        """
        Sets the initial condition for a vortex ring configuration in a 3-dimensional system
        """
        if not (self.dim == 3):
            raise Exception("The dimension of the system must be 3 for a vortex ring configuration.")

        if position is None:
            position = self.rmid

        if radius is None:
            radius = self.xmax / 3

        theta = self.calc_angle_field_vortex_ring(position=position, radius=radius, normal_vector=normal_vector)

        if self.psi is None:
            self.psi = 1

        self.psi = self.psi * np.exp(1j * theta)
        self.psi_f = np.fft.fftn(self.psi)

    def conf_vortex_remover(self, nodes, Area):
        '''
        Function that finds and removes vortices outside of the area defined by the corners
        (x1,y1), (x1,y2), (x2,y1), (x2,y2)
        Args:
            nodes (list) a list containing the vortices
            Area  (array) list on the format (x1,x2,y1,y2)
        returns:
        '''
        for vortex in nodes:
            x_coord = vortex['position'][0]
            y_coord = vortex['position'][1]
            if not (Area[0] < x_coord and x_coord < Area[1] \
                    and Area[2] < y_coord and y_coord < Area[3]):
                self.conf_insert_vortex(charge=-1 * vortex['charge'], position=[x_coord + self.dx, y_coord])
                # self.conf_insert_vortex(charge=vortex['charge'], position=[7, 0])

    def conf_dissipative_frame(self, d=7, wx=50, wy=50, wz=50):
        '''
        This function sets self.gamma so that it has a low value in the bulk and a large value near the edges.
        This sets a dissipative frame around the computational domain
         Args
             d (float): length of the interface between the low gamma and high gamma regions
            wx (float): distance fom center to the frame in x-direction
            wy (float):    -- " --                         y-direction
            wz (float):     -- " --                         z-direction
        return:
            modify self.gamma
        '''
        if self.dim == 2:
            X, Y = np.meshgrid(self.x, self.y, indexing='ij')
            gammax = self.gamma + 1 / 2 * (2 + np.tanh((X - self.xmid - wx) / d) - np.tanh((X - self.xmid + wx) / d))
            gammay = self.gamma + 1 / 2 * (2 + np.tanh((Y - self.ymid - wy) / d) - np.tanh((Y - self.ymid + wy) / d))
            self.gamma = np.real(np.maximum(gammax, gammay))
        elif self.dim == 3:
            X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')
            gammax = self.gamma + 1 / 2 * (2 + np.tanh((X - self.xmid - wx) / d) - np.tanh((X - self.xmid + wx) / d))
            gammay = self.gamma + 1 / 2 * (2 + np.tanh((Y - self.ymid - wy) / d) - np.tanh((Y - self.ymid + wy) / d))
            gammaz = self.gamma + 1 / 2 * (2 + np.tanh((Z - self.zmid - wz) / d) - np.tanh((Z - self.zmid + wz) / d))
            self.gamma = np.real(np.maximum(gammax, gammay, gammaz))
        else:
            raise Exception("This feature is not yet available for the given dimension.")


    # Time evolution
    def evolve_dGPE(self, number_of_steps, method='ETD2RK'):
        '''
       Evolver for the dGPE.
           Args:
               number_of_steps (int) the number of time steps that we are evolving the equation
               method (string, optional) the integration method we want to use. ETD2RK is sett as default
           returns:
               Updates the self.psi and self.psi_f
       '''

        k2 = self.calc_k2()
        omega_f = (1j + self.gamma) * (1 - 1 / 2 * k2)

        integrating_factors_f, solver = self.calc_integrating_factors_f_and_solver(omega_f, method)

        for n in tqdm(range(number_of_steps), desc='evolving the dGPE'):
            self.psi, self.psi_f = solver(integrating_factors_f,
                                          self.calc_nonlinear_evolution_function_f,
                                          self.psi, self.psi_f)



    def evolve_relax(self, number_of_steps, method='ETD2RK'):
        '''
        Evolver for the dGPE in imaginary time that relax the equation closer to the ground state
            Args:
                number_of_steps (int) the number of time steps that we are evolving the equation
                method (string, optional) the integration method we want to use. ETD2RK is sett as default
            returns:
                Updates the self.psi and self.psi_f
        '''
        temp_t = self.time
        gamma0 = self.gamma

        self.gamma = 1 - 1j

        print("Relaxing the BoseEinsteinCondensate...")
        self.evolve_dGPE(number_of_steps, method)

        self.gamma = gamma0
        self.time = temp_t

    def evolve_comoving_dGPE(self, number_of_steps, velx, method='ETD2RK'):
        '''
        Evolver for the dGPE in the comoving frame.
        This evolver assume that the stirring is in the x-direction and that gamma is spatialy dependent
            Args:
                number_of_steps (int) the number of time steps that we are evolving the equation
                velx (float) velocity in x direction
                method (string, optional) the integration method we want to use. ETD2RK is sett as default
            returns:
                Updates the fields self.psi and self.psi_f
                '''
        k2 = self.calc_k2()

        omega_f = (1j) * (1 - 1 / 2 * k2) + velx * self.dif[0]

        if method == 'ETD2RK':
            integrating_factors_f = self.calc_evolution_integrating_factors_ETD2RK(omega_f)
            solver = self.evolve_ETD2RK_loop
        elif method == 'ETD4RK':
            integrating_factors_f = self.calc_evolution_integrating_factors_ETD4RK(omega_f)
            solver = self.evolve_ETD4RK_loop
        else:
            raise Exception('This method is not implemented')

        for n in tqdm(range(number_of_steps), desc='evolving the dGPE in the comoving frame'):
            self.psi, self.psi_f = solver(integrating_factors_f, self.calc_nonlinear_evolution_term_comoving_f,
                                          self.psi, self.psi_f)

    # CALCULATION FUNCTIONS

    def calc_nonlinear_evolution_function_f(self, psi):
        """
        Calculates the non-linear evolution term of the dGPE
        Args:
             psi (numpy.ndarray): the wavefunction at a given time.
        returns:
            (numpy.ndarray): the non-linear evolution term
        """
        psi2 = np.abs(psi) ** 2
        return np.fft.fftn((1j + self.gamma) * (-self.V_ext() - psi2) * psi)

    def calc_nonlinear_evolution_term_comoving_f(self, psi):
        """
               Calculates the non-linear evolution term of the dGPE when gamma is not a constant.
               Relevant for example in the comoving frame when we have a dissipative frame around the edge.
               Args:
                   psi (numpy.ndarray): the wavefunction at a given time.
               Returns:
                    (numpy.ndarray): the non-linear evolution term
               """
        psi2 = np.abs(psi) ** 2
        term1 = np.fft.fftn(-(1j + self.gamma) * (self.V_ext() + psi2) * psi)
        term2 = np.fft.fftn(self.gamma * psi)
        term3 = 0.5 * np.fft.fftn(self.gamma * np.fft.ifftn(-self.calc_k2() * np.fft.fftn(psi)))
        return (term1 + term2 + term3)

        # Functions for callculationg properties of the BoseEinsteinCondensate

    def calc_superfluid_current(self):
        """
        Function that calculates the superfluid current
        retuns:
            (numpy.ndarray) the superfluid current
        """
        if self.dim == 2:
            J_s = np.zeros((self.dim,self.xRes,self.yRes))
        elif self.dim == 3:
            J_s = np.zeros((self.dim, self.xRes, self.yRes,self.zRes))
        else:
            raise(Exception('Calculation of the  superfluid current is not implemented in this dimension'))
        for i in range(self.dim):
            J_s[i] = np.imag( np.conj(self.psi) * np.fft.ifftn(1j*self.k[i] *self.psi_f ))
        return J_s

    def calc_velocity(self):
        """
        calculates the weighted velocity field
        returns:
            (numpy.ndarray) the weighted velocity field
        """
        if self.dim == 2:
            u = np.zeros((self.dim, self.xRes, self.yRes))
        elif self.dim == 3:
            u = np.zeros((self.dim, self.xRes, self.yRes, self.zRes))
        else:
            raise (Exception('Calculation of the weighted velocity is not implemented in this dimension'))
        theta = np.angle(self.psi)
        for i in range(self.dim):
            u[i] = np.imag(np.exp(-1j*theta)* np.fft.ifftn(1j * self.k[i] * self.psi_f))
        return u

    def calc_kinetic_energy(self):
        """
        Calculates the kinetic energy.
        returns:
            (float) the kinetic energy
        """
        u = self.calc_velocity()
        u2 = np.sum(u[i]**2 for i in range(self.dim))
        return 0.5*self.calc_integrate_field(u2)


    def calc_hamiltonian_density(self):
        """
        Function that calculates the hamiltonian density
        returns:
            (numpy.ndarray) the hamiltonian density
        """
        k2 = self.calc_k2()
        interaction_term = 1/2*np.abs(self.psi)**4
        potential_term = (self.V_ext() - 1 )* np.abs(self.psi)**2
        laplacian_term = -1/2 *np.real( np.conj(self.psi) * np.fft.ifftn(-k2*self.psi_f))
        return laplacian_term +potential_term + interaction_term

    def calc_hamiltonian(self):
        """
        Function that calculates the Hamitlonian
        returns:
            (Float) the Hamiltonian
        """
        H = self.calc_hamiltonian_density()
        return self.calc_integrate_field(H)

    def calc_harmonic_potential(self, R_tf):
        """
        Set returns a harmonic trap with R_tf being the Thomas-Fermi radius
        Args:
                R_tf (float): The Thomas-Fermi radius
        returns:
                A harmonic potential
        """
        trapping_strength = 1 / (R_tf ** 2)
        if self.dim == 1:
            return trapping_strength * (self.x - self.xmid) ** 2
        if self.dim == 2:
            return trapping_strength * (((self.x - self.xmid) ** 2).reshape(self.xRes, 1)
                                        + ((self.y - self.ymid) ** 2).reshape(1, self.yRes))
        if self.dim == 3:
            return trapping_strength * (((self.x - self.xmid) ** 2).reshape(self.xRes, 1, 1)
                                        + ((self.y - self.ymid) ** 2).reshape(1, self.yRes, 1)
                                        + ((self.z - self.zmid) ** 2).reshape(1, 1, self.zRes))

    def calc_gaussian_stirring_potential(self, size, strength, position):
        """
        Function for calculate a gaussian potential
        Args:
            size (float) size of the stirrer
            strength (float) strength of the stirrer, i.e how large the potential is at the stirrers position
            position (array) the position of the stirrer
        returns:
            (numpy.ndarray) a gaussian potential
        """

        if self.dim == 1:
            return strength * np.exp(-(self.x - position[0]) ** 2 / size ** 2)

        elif self.dim == 2:
            return strength * np.exp(-(((self.x - position[0]) ** 2).reshape(self.xRes, 1)
                                       + ((self.y - position[1]) ** 2).reshape(1, self.yRes)) / (size ** 2))
        elif self.dim == 3:
            return strength * np.exp(-(((self.x - position[0]) ** 2).reshape(self.xRes, 1, 1)
                                       + ((self.y - position[1]) ** 2).reshape(1, self.yRes, 1)
                                       + ((self.z - position[2]) ** 2).reshape(1, 1, self.zRes)) / (size ** 2))

    def calc_force_on_external_potential(self):
        """ calculates the average force acting on the external potential.
        returns:
            (numpy.ndarray) average force on the potential
        """
        Force =np.zeros(self.dim)
        potential_f = np.fft.ifftn(self.V_ext())
        for i in range(self.dim):
            Force_density = np.real(np.abs(self.psi)**2 * np.fft.ifftn(1j*self.k[i]* potential_f))
            Force[i] = self.calc_integrate_field(Force_density)
        return Force
        #TODO: It is not clear to me exactly what this function does (Vidar 04.12.23)



### Functions for calculating vortex properties
    def calc_vortex_density(self, psi=None):

        if psi is None:
            psi = self.psi

        return self.calc_defect_density([np.real(psi), np.imag(psi)])

    def calc_vortex_density_singular(self):
        # TODO: Insert the correct value of the equilibrium of psi, based on theory (Vidar 03.12.23)
        return self.calc_defect_density([np.real(self.psi), np.imag(self.psi)])

    def calc_vortex_velocity_field(self, dt_psi, psi=None):

        if psi is None:
            psi = self.psi

        return self.calc_defect_velocity_field([np.real(psi), np.imag(psi)],
                                        [np.real(dt_psi), np.imag(dt_psi)])

    def calc_vortex_nodes(self, dt_psi=None):
        """
        Calculate the positions and charges of vortex nodes based on the defect density.
        Returns:
            list: A list of dictionaries representing the vortex nodes. Each dictionary contains the following keys:
                  - 'position_index': The position index of the vortex node in the defect density array.
                  - 'charge': The charge of the vortex node.
                  - 'position': The position of the vortex node as a list [x, y].
        """

        # Calculate defect density
        rho = self.calc_vortex_density(self.psi)

        if dt_psi is not None:

            velocity_field = self.calc_vortex_velocity_field(dt_psi, self.psi)

        vortex_nodes = []

        if self.dim == 2:
            # Parameters to tune to make the algorithm work
            charge_tolerance = 0.2

            # Calculate the point where defect density is largest
            rho_max_index = np.unravel_index(np.argmax(np.abs(rho)), rho.shape)

            # Integrate the defect density around this point (i.e. in a disk around)
            disk = self.calc_region_disk(position = [self.x.flatten()[rho_max_index[0]],self.y.flatten()[rho_max_index[1]]],
                                         radius=1)
            charge = self.calc_integrate_field(rho, disk)

            # self.plot_field(rho)
            # plt.show()

            X, Y = np.meshgrid(self.x, self.y, indexing='ij')

            while abs(charge) > charge_tolerance:
                vortex = {}
                vortex['position_index'] = rho_max_index
                vortex['charge'] = np.sign(charge) * np.ceil(np.abs(charge))
                x = np.sum(disk * abs(rho) * X) / np.sum(disk * abs(rho))
                y = np.sum(disk * abs(rho) * Y) / np.sum(disk * abs(rho))
                vortex['position'] = [x, y]
                if dt_psi is not None:
                    vortex['velocity'] = [velocity_field[0][rho_max_index], velocity_field[1][rho_max_index]]
                else:
                    vortex['velocity'] = [float('nan'), float('nan')]

                # Calculate the velocity

                vortex_nodes.append(vortex)

                rho[disk] = 0
                rho_max_index = np.unravel_index(np.argmax(np.abs(rho)), rho.shape)

                disk = self.calc_region_disk(position=[self.x.flatten()[rho_max_index[0]], self.y.flatten()[rho_max_index[1]]],
                                             radius=1)
                charge = self.calc_integrate_field(rho, disk)

        elif self.dim == 3:
            # Parameters to tune to make the algorithm work
            charge_tolerance = 0.5
            integration_radius = 2
            cylinder_height = 1

            rho_norm = np.sqrt(rho[0]**2+rho[1]**2+rho[2]**2)

            # Calculate the point where defect density is largest
            rho_max_index = np.unravel_index(np.argmax(rho_norm), rho_norm.shape)
            # Integrate the defect density around this point (i.e. in cylinder around)
            tangent_vector = np.array([rho[0][rho_max_index],rho[1][rho_max_index],rho[2][rho_max_index]])
            tangent_vector = tangent_vector/np.linalg.norm(tangent_vector)
            cylinder = self.calc_region_cylinder(position=[self.x.flatten()[rho_max_index[0]], self.y.flatten()[rho_max_index[1]], self.z.flatten()[rho_max_index[2]]],
                                                 radius = integration_radius,
                                                 normal_vector = tangent_vector,
                                                 height = cylinder_height)
            charge = self.calc_integrate_field(rho_norm, cylinder)/cylinder_height

            X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')

            while charge > charge_tolerance:
                vortex = {}
                vortex['position_index'] = rho_max_index
                vortex['tangent_vector'] = tangent_vector
                x = np.sum(cylinder * abs(rho_norm) * X) / np.sum(cylinder * abs(rho_norm))
                y = np.sum(cylinder * abs(rho_norm) * Y) / np.sum(cylinder * abs(rho_norm))
                z = np.sum(cylinder * abs(rho_norm) * Z) / np.sum(cylinder * abs(rho_norm))

                vortex['position'] = [x, y, z]

                if dt_psi is not None:
                    vortex['velocity'] = [velocity_field[0][rho_max_index],
                                          velocity_field[1][rho_max_index],
                                          velocity_field[2][rho_max_index]]
                else:
                    vortex['velocity'] = [float('nan'),
                                          float('nan'),
                                          float('nan')]
                vortex_nodes.append(vortex)

                rho_norm[cylinder] = 0

                # self.plot_field(rho_norm)
                # plt.draw()
                # plt.pause(0.05)

                rho_max_index = np.unravel_index(np.argmax(rho_norm), rho_norm.shape)
                tangent_vector = np.array([rho[0][rho_max_index], rho[1][rho_max_index], rho[2][rho_max_index]])
                tangent_vector = tangent_vector / np.linalg.norm(tangent_vector)

                cylinder = self.calc_region_cylinder(position=[self.x.flatten()[rho_max_index[0]], self.y.flatten()[rho_max_index[1]],
                                                               self.z.flatten()[rho_max_index[2]]],
                                                     radius=integration_radius,
                                                     normal_vector = tangent_vector,
                                                     height = cylinder_height)
                charge = self.calc_integrate_field(rho_norm, cylinder) / cylinder_height


        return vortex_nodes

    # Plot functions

    def plot_vortex_nodes(self, vortex_nodes, ax=None):

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

            for vortex in vortex_nodes:

                if vortex['charge'] > 0:
                    x_coords_pos.append(vortex['position'][0])
                    y_coords_pos.append(vortex['position'][1])
                    vx_coords_pos.append(vortex['velocity'][0])
                    vy_coords_pos.append(vortex['velocity'][1])
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
            ax.set_aspect('equal')
            ax.set_facecolor('none')

            ax.set_xlabel('$x/a_0$')
            ax.set_ylabel('$y/a_0$')

            ax.set_xlim([0, self.xmax-self.dx])
            ax.set_ylim([0, self.ymax-self.dy])

        elif self.dim == 3:
            # Plotting options
            quiver_scale = 2 # The scale of the quiver arrows

            if ax == None:
                ax = plt.gcf().add_subplot(111, projection='3d')
               # ax = plt.gca()

            x_coords = []
            y_coords = []
            z_coords = []

            tx = []
            ty = []
            tz = []

            vx = []
            vy = []
            vz = []


            for vortex in vortex_nodes:
                x_coords.append(vortex['position'][0])
                y_coords.append(vortex['position'][1])
                z_coords.append(vortex['position'][2])

                tx.append(vortex['tangent_vector'][0])
                ty.append(vortex['tangent_vector'][1])
                tz.append(vortex['tangent_vector'][2])

                vx.append(vortex['velocity'][0])
                vy.append(vortex['velocity'][1])
                vz.append(vortex['velocity'][2])

            tx = np.array(tx)
            ty = np.array(ty)
            tz = np.array(tz)

            vx = np.array(vx)
            vy = np.array(vy)
            vz = np.array(vz)

            if not len(vx) == 0:
                v2 =vx**2 + vy**2 + vz**2
                v_norm = np.sqrt(max(v2))
            else:
                v_norm = 1

            #ax.scatter(x_coords, y_coords, z_coords, marker='o', color='black')
            ax.quiver(x_coords, y_coords, z_coords, quiver_scale*tx, quiver_scale*ty, quiver_scale*tz, color='blue')
            ax.quiver(x_coords, y_coords, z_coords, quiver_scale*vx/v_norm, quiver_scale*vy/v_norm, quiver_scale*vz/v_norm, color='green')

            ax.set_xlabel('$x/a_0$')
            ax.set_ylabel('$y/a_0$')
            ax.set_zlabel('$z/a_0$')

            ax.set_xlim([0, self.xmax-self.dx])
            ax.set_ylim([0, self.ymax-self.dy])
            ax.set_zlim([0, self.zmax-self.dz])

            ax.set_aspect('equal')
        ax.grid(True)

    
