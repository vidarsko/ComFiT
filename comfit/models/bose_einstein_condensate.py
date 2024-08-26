from typing import Callable

import numpy as np
import matplotlib.pyplot as plt
from comfit.core.base_system import BaseSystem
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D  
import scipy as sp

class BoseEinsteinCondensate(BaseSystem):
    def __init__(self, dimension: int, **kwargs):
        """Initializes a system to simulate a Bose-Einstein Condensate using the Gross-Pitaevskii equation.

        Args:
            dimension : int
                The dimension of the system.
            kwargs : dict, optional
                Optional keyword arguments to set additional parameters. see
                https://comfitlib.com/ClassBoseEinsteinCondensate/

        Returns:
            The system object representing the BoseEinsteinCondensate simulation. BoseEinsteinCondensate object

        Example:
            bec = BoseEinsteinCondensate(3,xRes=101,yRes=101,zRes=101, gamma=0.5)
            Creates a BoseEinsteinCondensate system with 3 dimensions and a spatial resolution of 101 in all directions.
            The dissipative factor gamma is set to 0.5.
        """

        # First initialize the BaseSystem class
        super().__init__(dimension, **kwargs)

        # Type of the system
        self.psi = None
        self.psi_f = None
        self.type = 'BoseEinsteinCondensate'

        # Default simulation parameters
        self.gamma = 0.01 if 'gamma' not in kwargs else kwargs['gamma']  # Dissipation (gamma)

        self.V_ext = lambda t: 0  # External potential, note this is a function of t

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __str__(self) -> str:
        """Output a string representation of the system.

        Args:
            None

        Returns:
            A string representation of the system.
        """
        return f"ComFiT object: BoseEinsteinCondensate\n \
                Dimension: {self.dim}\n"

    # CONFIGURATION FUNCTIONS
    def conf_initial_condition_disordered(self, noise_strength: float = 0.01) -> np.ndarray:
        """Sets disordered initial condition for the BoseEinsteinCondensate with some thermal flcutiations
        
        Args:
            noise_strength: the strength of the noise

        Returns:
            Sets the value of self.psi and self.psi_f
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
        self.psi_f = sp.fft.fftn(self.psi)

    def conf_external_potential(self, V_ext: Callable | float, additive: bool = False) -> None:
        """Sets the external potential of the system.

        Args:
            V_ext (function or float): the external potential
            additive (bool, optional): whether to add the new potential to the existing potential or not

        Returns:
            Modifies the value of self.V_ext
        """

        if not callable(V_ext):
            original_V_ext = V_ext  # Preserve the original value of V_ext
            V_ext = lambda t: original_V_ext
        
        if additive:
            self.V_ext = lambda t: self.V_ext(t) + V_ext(t)
        else:
            self.V_ext = V_ext



    def conf_initial_condition_Thomas_Fermi(self) -> None:
        """Finds the Thomas_Fermi ground state.

        Must be precided by an energy relaxation to find the true ground state

        Args:
            None
        
        Returns: 
            Sets the value of self.psi and self.psi_f
        """
        V_0 = np.zeros(self.dims) + self.V_ext(self.time)
        self.psi = np.emath.sqrt(1 - V_0)

        self.psi[V_0 > 1] = 0
        self.psi_f = sp.fft.fftn(self.psi)

    # CONFIGURATION FUNCTIONS
    def conf_insert_vortex(self, charge: int = 1, position: list[float] | None = None):
        """Sets the initial condition for a vortex dipole

        Args:
            charge (int): the charge of the vortex
            position (list): the position of the vortex

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
        self.psi_f = sp.fft.fftn(self.psi)

    def conf_insert_vortex_dipole(
        self,
        dipole_vector: list[float] | None = None,
        dipole_position: list[float] | None = None
    ):
        """Sets the initial condition for a vortex dipole configuration in a 2-dimensional system.

        Args:
            None

        Returns:
            Modifies the value of self.psi and self.psi_f

        Raises:
            Exception: If the dimension of the system is not 2.
        """
        if not (self.dim == 2):
            raise Exception("The dimension of the system must be 2 for a vortex dipole configuration.")

        if self.psi is None:
            self.conf_initial_condition_Thomas_Fermi()

        if dipole_vector is None:
            dipole_vector = [(self.xmax-self.xmin) / 3, 0]

        if dipole_position is None:
            dipole_position = self.rmid

        self.psi = self.psi * np.exp(1j * self.calc_angle_field_vortex_dipole(dipole_vector, dipole_position))
        self.psi_f = sp.fft.fftn(self.psi)

    def conf_insert_vortex_ring(
        self,
        position: list[float] | None = None,
        radius: float | None = None,
        normal_vector: list[float] = [0, 0, 1]
    ) -> None:
        """Sets the initial condition for a vortex ring configuration in a 3-dimensional system
        
        Args:
            position: the position of the vortex ring (list)
            radius: the radius of the vortex ring (float)
            normal_vector: the normal vector of the vortex ring (list)
        
        Returns:
            Modifies the value of self.psi and self.psi_f
        """
        if not (self.dim == 3):
            raise Exception("The dimension of the system must be 3 for a vortex ring configuration.")

        if position is None:
            position = self.rmid

        if radius is None:
            radius = (self.xmax -self.xmin)/ 3

        theta = self.calc_angle_field_vortex_ring(position=position, radius=radius, normal_vector=normal_vector)

        if self.psi is None:
            self.conf_initial_condition_Thomas_Fermi()

        self.psi = self.psi * np.exp(1j * theta)
        self.psi_f = sp.fft.fftn(self.psi)

    def conf_vortex_remover(self, nodes: list[dict], Area: list[float]) -> None:
        '''Removes vortices
        
        Function that finds and removes vortices outside of the area defined by the corners
        (x1,y1), (x1,y2), (x2,y1), (x2,y2)

        Args:
            nodes (list) a list containing the vortices
            Area  (array) list on the format (x1,x2,y1,y2)

        Returns:
            Modifies the value of self.psi and self.psi_f
        '''
        for vortex in nodes:
            x_coord = vortex['position'][0]
            y_coord = vortex['position'][1]
            if not (Area[0] < x_coord and x_coord < Area[1] \
                    and Area[2] < y_coord and y_coord < Area[3]):
                self.conf_insert_vortex(charge=-1 * vortex['charge'], position=[x_coord + self.dx, y_coord])
                # self.conf_insert_vortex(charge=vortex['charge'], position=[7, 0])

    def conf_dissipative_frame(self, d: float = 7, wx: float = 50, wy: float = 50, wz: float = 50):
        '''Configures a dissipative frame around the computational domain

        This function sets self.gamma so that it has a low value in the bulk and a large value near the edges.
        This sets a dissipative frame around the computational domain
        
        Args:
            d: length of the interface between the low gamma and high gamma regions (float)
            wx: distance fom center to the frame in x-direction (float)
            wy:    -- " --                         y-direction (float)
            wz:     -- " --                         z-direction (float)
        
        Returns:
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
    def evolve_dGPE(self, number_of_steps: int , method: str = 'ETD2RK') -> None:
        '''Evolver for the dGPE.
        
        Args:
            number_of_steps: the number of time steps that we are evolving the equation  (int)
            method: the integration method we want to use. ETD2RK is sett as default  (string, optional)
        
        Returns:
            Updates the self.psi and self.psi_f
       '''

        k2 = self.calc_k2()
        omega_f = (1j + self.gamma) * (1 - 1 / 2 * k2)

        integrating_factors_f, solver = self.calc_integrating_factors_f_and_solver(omega_f, method)

        for n in tqdm(range(number_of_steps), desc='evolving the dGPE'):
            self.psi, self.psi_f = solver(integrating_factors_f,
                                          self.calc_nonlinear_evolution_function_f,
                                          self.psi, self.psi_f)

    def evolve_relax(self, number_of_steps: int, method: str = 'ETD2RK') -> None:
        '''Evolver for the dGPE in imaginary time that relax the equation closer to the ground state
        
        Args:
            number_of_steps: the number of time steps that we are evolving the equation (int)
            method: the integration method we want to use. ETD2RK is sett as default (string)
        
        Returns:
            Updates the self.psi and self.psi_f
        '''
        temp_t = self.time
        gamma0 = self.gamma

        self.gamma = 1 - 1j

        temp_V = self.V_ext

        self.conf_external_potential(temp_V(temp_t))

        print("Relaxing the BoseEinsteinCondensate...")
        self.evolve_dGPE(number_of_steps, method)

        self.gamma = gamma0
        self.time = temp_t
        
        self.conf_external_potential(temp_V)

    def evolve_comoving_dGPE(self, number_of_steps: int, velx: float, method: str = 'ETD2RK') -> None:
        '''Evolver for the dGPE in the comoving frame.

        This evolver assume that the stirring is in the x-direction and that gamma is spatialy dependent
        
        Args:
            number_of_steps: the number of time steps that we are evolving the equation (int)
            velx: velocity in x direction (float) 
            method: the integration method we want to use. ETD2RK is sett as default (string)
        
        Returns:
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

    def calc_nonlinear_evolution_function_f(self, psi: np.ndarray, t) -> np.ndarray:
        """Calculates the non-linear evolution term of the dGPE
        
        Args:
            psi: the wavefunction at a given time. (numpy.ndarray)
        
        Returns:
            The non-linear evolution term (numpy.ndarray)
        """
        
        psi2 = np.abs(psi) ** 2
        return sp.fft.fftn((1j + self.gamma) * (-self.V_ext(t) - psi2) * psi)

    def calc_nonlinear_evolution_term_comoving_f(self, psi: np.ndarray, t) -> np.ndarray:
        """Calculates the non-linear evolution term of the dGPE when gamma is not a constant.

        Relevant for example in the comoving frame when we have a dissipative frame around the edge.
        
        Args:
            psi: the wavefunction at a given time. (numpy.ndarray)
        
        Returns:
            the non-linear evolution term (numpy.ndarray)
        """
        psi2 = np.abs(psi) ** 2
        term1 = sp.fft.fftn(-(1j + self.gamma) * (self.V_ext(t) + psi2) * psi)
        term2 = sp.fft.fftn(self.gamma * psi)
        term3 = 0.5 * sp.fft.fftn(self.gamma * sp.fft.ifftn(-self.calc_k2() * sp.fft.fftn(psi)))
        return (term1 + term2 + term3)

        # Functions for callculationg properties of the BoseEinsteinCondensate

    def calc_superfluid_current(self) -> np.ndarray:
        """Calculates the superfluid current

        Args:
            None
        
        Returns:
            The superfluid current (numpy.ndarray) 
        """
        if self.dim == 2:
            J_s = np.zeros((self.dim,self.xRes,self.yRes))
        elif self.dim == 3:
            J_s = np.zeros((self.dim, self.xRes, self.yRes,self.zRes))
        else:
            raise(Exception('Calculation of the  superfluid current is not implemented in this dimension'))
        for i in range(self.dim):
            J_s[i] = np.imag( np.conj(self.psi) * sp.fft.ifftn(1j*self.k[i] *self.psi_f ))
        return J_s

    def calc_velocity(self) -> np.ndarray:
        """Calculates the weighted velocity field
        
        Args:
            None

        Returns:
            The weighted velocity field (numpy.ndarray) 
        """
        if self.dim == 2:
            u = np.zeros((self.dim, self.xRes, self.yRes))
        elif self.dim == 3:
            u = np.zeros((self.dim, self.xRes, self.yRes, self.zRes))
        else:
            raise (Exception('Calculation of the weighted velocity is not implemented in this dimension'))
        theta = np.angle(self.psi)
        for i in range(self.dim):
            u[i] = np.imag(np.exp(-1j*theta)* sp.fft.ifftn(1j * self.k[i] * self.psi_f))
        return u

    def calc_kinetic_energy(self) -> float:
        """Calculates the kinetic energy.

        Args:
            None

        Returns:
            The kinetic energy (float) 
        """
        u = self.calc_velocity()
        u2 = np.sum(u[i]**2 for i in range(self.dim))
        return 0.5*self.calc_integrate_field(u2)


    def calc_hamiltonian_density(self) -> np.ndarray:
        """Calculates the hamiltonian density

        Args:
            None

        Returns:
            The hamiltonian density (numpy.ndarray) 
        """
        k2 = self.calc_k2()
        interaction_term = 1/2*np.abs(self.psi)**4
        potential_term = (self.V_ext(self.time) - 1 )* np.abs(self.psi)**2
        laplacian_term = -1/2 *np.real( np.conj(self.psi) * sp.fft.ifftn(-k2*self.psi_f))
        return laplacian_term +potential_term + interaction_term

    def calc_hamiltonian(self) -> float:
        """Function that calculates the Hamiltonian

        Args:
            None

        Returns:
            The Hamiltonian
        """
        H = self.calc_hamiltonian_density()
        return self.calc_integrate_field(H)

    def calc_harmonic_potential(self, R_tf: float) -> np.ndarray:
        """Calculates a harmonic trap with R_tf being the Thomas-Fermi radius
        
        Args:
            R_tf (float): The Thomas-Fermi radius
        
        Returns:
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

    def calc_gaussian_stirring_potential(self, size: float, strength: float, position: np.ndarray) -> np.ndarray:
        """Calculates a gaussian potential
        
        Args:
            size (float) size of the stirrer
            strength (float) strength of the stirrer, i.e how large the potential is at the stirrers position
            position (array) the position of the stirrer
       
        Retuns:
            A gaussian potential (numpy.ndarray) 
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

    def calc_force_on_external_potential(self) -> np.ndarray: 
        """Calculates the average force acting on the external potential.
        
        Args:
            None
        
        Returns:
            Average force on the potential (numpy.ndarray) 
        """
        Force =np.zeros(self.dim)
        potential_f = sp.fft.ifftn(self.V_ext(self.time))
        for i in range(self.dim):
            Force_density = np.real(np.abs(self.psi)**2 * sp.fft.ifftn(1j*self.k[i]* potential_f))
            Force[i] = self.calc_integrate_field(Force_density)
        return Force
        #TODO: It is not clear to me exactly what this function does (Vidar 04.12.23)



    ## Functions for calculating vortex properties
    def calc_vortex_density(self, psi: np.ndarray | None = None) -> np.ndarray:
        """Calculates the vortex density of the system.

        Args:
            psi: The wavefunction of the system. (numpy.ndarray)
        
        Returns:
            The vortex density of the system. numpy.ndarray
        """
        if psi is None:
            psi = self.psi

        return self.calc_defect_density([np.real(psi), np.imag(psi)])

    def calc_vortex_density_singular(self) -> np.ndarray:
        """Calculates the vortex density of the system using the singular method.

        Args:
            None
        
        Returns:
            The vortex density of the system. (numpy.ndarray)
        """
        # TODO: Insert the correct value of the equilibrium of psi, based on theory (Vidar 03.12.23)
        return self.calc_defect_density([np.real(self.psi), np.imag(self.psi)])

    def calc_vortex_velocity_field(self, dt_psi: np.ndarray, psi: np.ndarray | None = None) -> np.ndarray:
        """Calculates the vortex velocity field of the system.

        Args:
            dt_psi: The time derivative of the wavefunction of the system. (numpy.ndarray)
            psi: The wavefunction of the system. (numpy.ndarray)

        Returns:
            The vortex velocity field of the system. numpy.ndarray: 
        """
        if psi is None:
            psi = self.psi

        return self.calc_defect_velocity_field([np.real(psi), np.imag(psi)],
                                        [np.real(dt_psi), np.imag(dt_psi)])

    def calc_vortex_nodes(self, dt_psi: np.ndarray | None = None) -> list[dict]:
        """
        Calculate the positions and charges of vortex nodes based on the defect density.

        Args:
            dt_psi: The time derivative of the wavefunction of the system.     (numpy.ndarray)
    
        Returns:
            List of dictionaries representing the vortex nodes. Each dictionary contains the following keys:
                  - 'position_index': The position index of the vortex node in the defect density array.
                  - 'charge': The charge of the vortex node.
                  - 'position': The position of the vortex node as a list [x, y].
                  - 'velocity': The velocity of the vortex node as a list [vx, vy].
        """

        # Calculate defect density
        rho = self.calc_vortex_density(self.psi)

        if dt_psi is not None:
            velocity_field = self.calc_vortex_velocity_field(dt_psi, self.psi)

        if self.dim == 2:
            vortex_nodes = self.calc_defect_nodes(np.abs(rho), 
                                                    charge_tolerance = 0.2,
                                                    integration_radius = self.a0)
            for vortex in vortex_nodes:
                vortex['charge'] = np.sign(rho[vortex['position_index']])
                if dt_psi is not None:
                    vortex['velocity'] = [velocity_field[0][vortex['position_index']], 
                                        velocity_field[1][vortex['position_index']]]
                else:
                    vortex['velocity'] = [float('nan'), float('nan')]
        elif self.dim == 3:
            vortex_nodes = self.calc_defect_nodes(np.sqrt(rho[0]**2 + rho[1]**2 + rho[2]**2), 
                                                charge_tolerance = 2*self.a0,
                                                integration_radius = 2*self.a0
                )
            for vortex in vortex_nodes:
                tangent_vector = np.array([rho[i][vortex['position_index']] for i in range(3)]),
                vortex['tangent_vector'] = tangent_vector[0]/np.linalg.norm(tangent_vector)
                
                if dt_psi is not None:
                    vortex['velocity'] = [velocity_field[0][vortex['position_index']], 
                                        velocity_field[1][vortex['position_index']], 
                                        velocity_field[2][vortex['position_index']]]
                else:
                    vortex['velocity'] = [float('nan'), float('nan'), float('nan')]

        return vortex_nodes

    # Plot functions

    def plot_vortex_nodes(self, vortex_nodes: list[dict], **kwargs):
        """Plots the vortex nodes in the system.

        Args:
            vortex_nodes (list): A list of dictionaries representing the vortex nodes. Each dictionary contains the following keys:
                                 - 'position_index': The position index of the vortex node in the defect density array.
                                 - 'charge': The charge of the vortex node.
                                 - 'position': The position of the vortex node as a list [x, y].
                                 - 'velocity': The velocity of the vortex node as a list [vx, vy].
            -**kwargs: Keyword arguments for the plot.
                See https://comfitlib.com/ClassBaseSystem/
                for a full list of keyword arguments.
        
        Returns:
            The axes on which the vortex nodes are plotted. (matplotlib.axes.Axes: )
        """
        # Check if an axis object is provided
        fig = kwargs.get('fig', plt.gcf())
        ax = kwargs.get('ax', None)
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

        return ax

    
