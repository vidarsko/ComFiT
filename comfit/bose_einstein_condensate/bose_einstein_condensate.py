from typing import Callable, Union, Optional

import numpy as np
import matplotlib.pyplot as plt
from comfit.core.base_system import BaseSystem
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D  
import scipy as sp

class BoseEinsteinCondensate(BaseSystem):
    def __init__(self, dim: int, **kwargs):
        """Initializes a system to simulate a Bose-Einstein Condensate using the Gross-Pitaevskii equation.

        Parameters
        ----------
        dim : int
            The dimension of the system.
        kwargs : dict, optional
            Optional keyword arguments to set additional parameters. see
            https://comfitlib.com/ClassBoseEinsteinCondensate/

        Returns
        -------
        BoseEinsteinCondensate
            The system object representing the BoseEinsteinCondensate simulation.

        Examples
        --------
        >>> bec = BoseEinsteinCondensate(3,xRes=101,yRes=101,zRes=101, gamma=0.5)
        Creates a BoseEinsteinCondensate system with 3 dimensions and a spatial resolution of 101 in all directions.
        The dissipative factor gamma is set to 0.5.
        """

        # First initialize the BaseSystem class
        super().__init__(dim, **kwargs)

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

        Returns
        -------
        str
            A string representation of the system.
        """
        return f"ComFiT object: BoseEinsteinCondensate\n \
                Dimension: {self.dim}\n"

    # CONFIGURATION FUNCTIONS
    def conf_initial_condition_disordered(self, noise_strength: float = 0.01) -> np.ndarray:
        """Sets disordered initial condition for the BoseEinsteinCondensate with some thermal fluctuations

        Parameters
        ----------
        noise_strength : float
            the strength of the noise

        Returns
        -------
        None
            Sets the value of self.psi and self.psi_f. self.psi always
            carries a leading component axis, so self.psi[0] is set to the
            newly configured disordered state.
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

        self.psi = np.array([noise_strength * self.psi])
        self.psi_f = sp.fft.fftn(self.psi)

    def conf_external_potential(self, V_ext: Union[Callable, float], additive: bool = False) -> None:
        """Sets the external potential of the system.

        Parameters
        ----------
        V_ext : function or float
            the external potential
        additive : bool, optional
            whether to add the new potential to the existing potential or not

        Returns
        -------
        None
            Modifies the value of self.V_ext
        """

        if not callable(V_ext):
            original_V_ext = V_ext  # Preserve the original value of V_ext
            V_ext = lambda t: original_V_ext
        
        if additive:
            self.V_ext = lambda t: self.V_ext(t) + V_ext(t)
        else:
            self.V_ext = V_ext



    def conf_initial_condition_thomas_fermi(self) -> None:
        """Finds the Thomas_Fermi ground state.

        Must be preceded by an energy relaxation to find the true ground state

        Returns
        -------
        None
            Sets the value of self.psi and self.psi_f. self.psi always
            carries a leading component axis, so self.psi[0] is set to the
            newly configured Thomas-Fermi state.
        """
        V_0 = np.zeros(self.dims) + self.V_ext(self.time)
        psi = np.emath.sqrt(1 - V_0)

        psi[V_0 > 1] = 0
        self.psi = np.array([psi])
        self.psi_f = sp.fft.fftn(self.psi)

    # CONFIGURATION FUNCTIONS
    def conf_insert_vortex(self, charge: int = 1, position: Optional[list[float]] = None):
        """Sets the initial condition for a vortex dipole

        Parameters
        ----------
        charge : int
            the charge of the vortex
        position : list
            the position of the vortex

        Returns
        -------
        None
            Modifies the value of self.psi and self.psi_f
        """
        if not (self.dim == 2):
            raise Exception("The dimension of the system must be 2 for a single point vortex.")

        if position is None:
            position = [self.xmid, self.ymid]

        if self.psi is None:
            self.conf_initial_condition_thomas_fermi()

        self.psi = self.psi * np.exp(1j * self.calc_angle_field_single_vortex(position=position, charge=charge))
        self.psi_f = sp.fft.fftn(self.psi)

    def conf_insert_vortex_dipole(
        self,
        dipole_vector: Optional[list[float]] = None,
        dipole_position: Optional[list[float]] = None
    ):
        """Sets the initial condition for a vortex dipole configuration in a 2-dimensional system.

        Returns
        -------
        None
            Modifies the value of self.psi and self.psi_f

        Raises
        ------
        Exception
            If the dimension of the system is not 2.
        """
        if not (self.dim == 2):
            raise Exception("The dimension of the system must be 2 for a vortex dipole configuration.")

        if self.psi is None:
            self.conf_initial_condition_thomas_fermi()

        if dipole_vector is None:
            dipole_vector = [self.size_x / 3, 0]

        if dipole_position is None:
            dipole_position = self.rmid

        self.psi = self.psi * np.exp(1j * self.calc_angle_field_vortex_dipole(dipole_vector, dipole_position))
        self.psi_f = sp.fft.fftn(self.psi)

    def conf_insert_vortex_filament(
            self,
            position = None,
            charge = None
    ) -> None:
        """ Insert a vortex line into the condensate. The vortex line is assumed to be elongated along the z-axis.

        Parameters
        ----------
        position : array_like, optional
            the position in the xy-plane of the vortex filament (vector)
        charge : int, optional
            the charge of the vortex filament

        Returns
        -------
        None
            Modifies the value of self.psi and self.psi_f
        """

        if not (self.dim == 3):
            raise Exception("The dimension of the system must be 3 for a vortex line configuration.")

        if self.psi is None:
            self.conf_initial_condition_thomas_fermi()
        
        if position is None:
            position = [self.xmid , self.ymid]

        if charge is None:
            charge = 1

        theta = self.calc_angle_field_vortex_filament(charge=charge,position=position)

        self.psi = self.psi*np.exp(1j *theta)
        self.psi_f = sp.fft.fftn(self.psi)

    def conf_insert_vortex_ring(
        self,
        position: Optional[list[float]] = None,
        radius: Optional[float] = None,
        normal_vector: list[float] = [0, 0, 1]
    ) -> None:
        """Sets the initial condition for a vortex ring configuration in a 3-dimensional system

        Parameters
        ----------
        position : list, optional
            the position of the vortex ring
        radius : float, optional
            the radius of the vortex ring
        normal_vector : list, optional
            the normal vector of the vortex ring

        Returns
        -------
        None
            Modifies the value of self.psi and self.psi_f
        """
        if not (self.dim == 3):
            raise Exception("The dimension of the system must be 3 for a vortex ring configuration.")

        if position is None:
            position = self.rmid

        if radius is None:
            radius = self.size_x/ 3

        theta = self.calc_angle_field_vortex_ring(position=position, radius=radius, normal_vector=normal_vector)

        if self.psi is None:
            self.conf_initial_condition_thomas_fermi()

        self.psi = self.psi * np.exp(1j * theta)
        self.psi_f = sp.fft.fftn(self.psi)

    def conf_vortex_remover(self, nodes: list[dict], area: list[float]) -> None:
        '''Removes vortices

        Function that finds and removes vortices outside of the area defined by the corners
        (x1,y1), (x1,y2), (x2,y1), (x2,y2)

        Parameters
        ----------
        nodes : list
            a list containing the vortices
        area : array
            list on the format (x1,x2,y1,y2)

        Returns
        -------
        None
            Modifies the value of self.psi and self.psi_f
        '''
        for vortex in nodes:
            x_coord = vortex['position'][0]
            y_coord = vortex['position'][1]
            if not (area[0] < x_coord and x_coord < area[1] \
                    and area[2] < y_coord and y_coord < area[3]):
                self.conf_insert_vortex(charge=-1 * vortex['charge'], position=[x_coord + self.dx, y_coord])
                # self.conf_insert_vortex(charge=vortex['charge'], position=[7, 0])

    def conf_dissipative_frame(self, interface_width: float = 7, 
                                    frame_width_x: float = None, 
                                    frame_width_y: float = None, 
                                    frame_width_z: float = None) -> None:

        '''Configures a dissipative frame around the computational domain

        This function sets self.gamma so that it has a low value in the bulk and a large value near the edges.
        This sets a dissipative frame around the computational domain

        Parameters
        ----------
        interface_width : float
            length of the interface between the low gamma and high gamma regions
        frame_width_x : float
            distance from center to the frame in x-direction
        frame_width_y : float
            -- " --                         y-direction
        frame_width_z : float
            -- " --                         z-direction

        Returns
        -------
        None
            modify self.gamma
        '''

        # If configuration is called for the first time, set initial_gamma
        if not hasattr(self, 'initial_gamma'):
            self.initial_gamma = self.gamma

        # Set default values for frame_width
        if frame_width_x is None:
            frame_width_x = 0.8*self.size_x
        if frame_width_y is None:
            frame_width_y = 0.8*self.size_y
        if frame_width_z is None:
            frame_width_z = 0.8*self.size_z

        # Set the frame
        if self.dim == 2:
            X, Y = np.meshgrid(self.x, self.y, indexing='ij')
            gammax = self.initial_gamma + 0.5 * (2 + np.tanh((X - self.xmid - frame_width_x/2) /interface_width) - np.tanh((X - self.xmid + frame_width_x/2) /interface_width))
            gammay = self.initial_gamma + 0.5 * (2 + np.tanh((Y - self.ymid - frame_width_y/2) /interface_width) - np.tanh((Y - self.ymid + frame_width_y/2) /interface_width))
            self.gamma = np.real(np.maximum(gammax, gammay))
        elif self.dim == 3:
            X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')
            gammax = self.initial_gamma + 0.5 * (2 + np.tanh((X - self.xmid - frame_width_x/2) /interface_width) - np.tanh((X - self.xmid + frame_width_x/2) /interface_width))
            gammay = self.initial_gamma + 0.5 * (2 + np.tanh((Y - self.ymid - frame_width_y/2) /interface_width) - np.tanh((Y - self.ymid + frame_width_y/2) /interface_width))
            gammaz = self.initial_gamma + 0.5 * (2 + np.tanh((Z - self.zmid - frame_width_z/2) /interface_width) - np.tanh((Z - self.zmid + frame_width_z/2) /interface_width))
            self.gamma = np.real(np.maximum(gammax, gammay, gammaz))
        else:
            raise Exception("This feature is not yet available for the given dimension.")


    # Time evolution
    def evolve_dGPE(self, number_of_steps: int , method: str = 'ETD2RK') -> None:
        '''Evolver for the dGPE.

        Parameters
        ----------
        number_of_steps : int
            the number of time steps that we are evolving the equation
        method : string, optional
            the integration method we want to use. ETD2RK is sett as default

        Returns
        -------
        None
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

        Parameters
        ----------
        number_of_steps : int
            the number of time steps that we are evolving the equation
        method : string
            the integration method we want to use. ETD2RK is sett as default

        Returns
        -------
        None
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

    def evolve_comoving_dGPE(self, number_of_steps: int, vel_x: float, method: str = 'ETD2RK') -> None:
        '''Evolver for the dGPE in the comoving frame.

        This evolver assume that the stirring is in the x-direction and that gamma is spatialy dependent

        Parameters
        ----------
        number_of_steps : int
            the number of time steps that we are evolving the equation
        vel_x : float
            velocity in x direction
        method : string
            the integration method we want to use. ETD2RK is sett as default

        Returns
        -------
        None
            Updates the fields self.psi and self.psi_f
        '''
        k2 = self.calc_k2()

        omega_f = (1j) * (1 - 1 / 2 * k2) + vel_x * self.dif[0]

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

        Parameters
        ----------
        psi : numpy.ndarray
            the wavefunction at a given time.

        Returns
        -------
        numpy.ndarray
            The non-linear evolution term
        """
        
        psi2 = np.abs(psi) ** 2
        return sp.fft.fftn((1j + self.gamma) * (-self.V_ext(t) - psi2) * psi)

    def calc_nonlinear_evolution_term_comoving_f(self, psi: np.ndarray, t) -> np.ndarray:
        """Calculates the non-linear evolution term of the dGPE when gamma is not a constant.

        Relevant for example in the comoving frame when we have a dissipative frame around the edge.

        Parameters
        ----------
        psi : numpy.ndarray
            the wavefunction at a given time.

        Returns
        -------
        numpy.ndarray
            the non-linear evolution term
        """
        psi2 = np.abs(psi) ** 2
        term1 = sp.fft.fftn(-(1j + self.gamma) * (self.V_ext(t) + psi2) * psi)
        term2 = sp.fft.fftn(self.gamma * psi)
        term3 = 0.5 * sp.fft.fftn(self.gamma * sp.fft.ifftn(-self.calc_k2() * sp.fft.fftn(psi)))
        return (term1 + term2 + term3)

        # Functions for callculationg properties of the BoseEinsteinCondensate

    def calc_superfluid_current(self) -> np.ndarray:
        """Calculates the superfluid current

        Returns
        -------
        numpy.ndarray
            The superfluid current
        """
        if self.dim == 2:
            J_s = np.zeros((self.dim,self.xRes,self.yRes))
        elif self.dim == 3:
            J_s = np.zeros((self.dim, self.xRes, self.yRes,self.zRes))
        else:
            raise(Exception('Calculation of the  superfluid current is not implemented in this dimension'))
        psi = self.psi[0]
        psi_f = self.psi_f[0]
        for i in range(self.dim):
            J_s[i] = np.imag( np.conj(psi) * sp.fft.ifftn(1j*self.k[i] *psi_f ))
        return J_s

    def calc_velocity(self) -> np.ndarray:
        """Calculates the weighted velocity field

        Returns
        -------
        numpy.ndarray
            The weighted velocity field
        """
        if self.dim == 2:
            u = np.zeros((self.dim, self.xRes, self.yRes))
        elif self.dim == 3:
            u = np.zeros((self.dim, self.xRes, self.yRes, self.zRes))
        else:
            raise (Exception('Calculation of the weighted velocity is not implemented in this dimension'))
        theta = np.angle(self.psi[0])
        psi_f = self.psi_f[0]
        for i in range(self.dim):
            u[i] = np.imag(np.exp(-1j*theta)* sp.fft.ifftn(1j * self.k[i] * psi_f))
        return u

    def calc_kinetic_energy(self) -> float:
        """Calculates the kinetic energy.

        Returns
        -------
        float
            The kinetic energy
        """
        u = self.calc_velocity()
        u2 = sum(u[i]**2 for i in range(self.dim))
        return 0.5*self.calc_integrate_field(u2)


    def calc_hamiltonian_density(self) -> np.ndarray:
        """Calculates the hamiltonian density

        Returns
        -------
        numpy.ndarray
            The hamiltonian density
        """
        psi = self.psi[0]
        psi_f = self.psi_f[0]
        k2 = self.calc_k2()
        interaction_term = 1/2*np.abs(psi)**4
        potential_term = (self.V_ext(self.time) - 1 )* np.abs(psi)**2
        laplacian_term = -1/2 *np.real( np.conj(psi) * sp.fft.ifftn(-k2*psi_f))
        return laplacian_term +potential_term + interaction_term

    def calc_hamiltonian(self) -> float:
        """Function that calculates the Hamiltonian

        Returns
        -------
        float
            The Hamiltonian
        """
        H = self.calc_hamiltonian_density()
        return self.calc_integrate_field(H)

    def calc_harmonic_potential(self, thomas_fermi_radius: float) -> np.ndarray:
        """Calculates a harmonic trap with thomas_fermi_radius being the Thomas-Fermi radius

        Parameters
        ----------
        thomas_fermi_radius : float
            The Thomas-Fermi radius

        Returns
        -------
        numpy.ndarray
            A harmonic potential
        """
        trapping_strength = 1 / (thomas_fermi_radius ** 2)
        if self.dim == 1:
            return trapping_strength * (self.x - self.xmid) ** 2
        if self.dim == 2:
            return trapping_strength * (((self.x - self.xmid) ** 2).reshape(self.xRes, 1)
                                        + ((self.y - self.ymid) ** 2).reshape(1, self.yRes))
        if self.dim == 3:
            return trapping_strength * (((self.x - self.xmid) ** 2).reshape(self.xRes, 1, 1)
                                        + ((self.y - self.ymid) ** 2).reshape(1, self.yRes, 1)
                                        + ((self.z - self.zmid) ** 2).reshape(1, 1, self.zRes))

    def calc_force_on_external_potential(self) -> np.ndarray: 
        """Calculates the average force acting on the external potential.

        Returns
        -------
        numpy.ndarray
            Average force on the potential
        """
        Force =np.zeros(self.dim)
        potential_f = sp.fft.ifftn(self.V_ext(self.time))
        for i in range(self.dim):
            Force_density = np.real(np.abs(self.psi[0])**2 * sp.fft.ifftn(1j*self.k[i]* potential_f))
            Force[i] = self.calc_integrate_field(Force_density)
        return Force
        #TODO: It is not clear to me exactly what this function does (Vidar 04.12.23)



    ## Functions for calculating vortex properties
    def calc_vortex_density(self, psi: Optional[np.ndarray] = None) -> np.ndarray:
        """Calculates the vortex density of the system.

        Parameters
        ----------
        psi : numpy.ndarray, optional
            The wavefunction of the system. If None, self.psi[0] is used.

        Returns
        -------
        numpy.ndarray
            The vortex density of the system.
        """
        if psi is None:
            psi = self.psi[0]

        return self.calc_defect_density([np.real(psi), np.imag(psi)])

    def calc_vortex_density_singular(self) -> np.ndarray:
        """Calculates the vortex density of the system using the singular method.

        Returns
        -------
        numpy.ndarray
            The vortex density of the system.
        """
        # TODO: Insert the correct value of the equilibrium of psi, based on theory (Vidar 03.12.23)
        return self.calc_defect_density([np.real(self.psi[0]), np.imag(self.psi[0])])

    def calc_vortex_velocity_field(self, dt_psi: np.ndarray, psi: Optional[np.ndarray] = None) -> np.ndarray:
        """Calculates the vortex velocity field of the system.

        Parameters
        ----------
        dt_psi : numpy.ndarray
            The time derivative of the wavefunction of the system.
        psi : numpy.ndarray, optional
            The wavefunction of the system. If None, self.psi[0] is used.

        Returns
        -------
        numpy.ndarray
            The vortex velocity field of the system.
        """
        if psi is None:
            psi = self.psi[0]

        return self.calc_defect_velocity_field([np.real(psi), np.imag(psi)],
                                        [np.real(dt_psi), np.imag(dt_psi)])

    def calc_vortex_nodes(self, dt_psi: Optional[np.ndarray] = None) -> list[dict]:
        """
        Calculate the positions and charges of vortex nodes based on the defect density.

        Parameters
        ----------
        dt_psi : numpy.ndarray, optional
            The time derivative of the wavefunction of the system.

        Returns
        -------
        list of dict
            List of dictionaries representing the vortex nodes. Each dictionary contains the following keys:
                  - 'position_index': The position index of the vortex node in the defect density array.
                  - 'charge': The charge of the vortex node.
                  - 'position': The position of the vortex node as a list [x, y].
                  - 'velocity': The velocity of the vortex node as a list [vx, vy].
        """

        # Calculate defect density
        rho = self.calc_vortex_density(self.psi[0])

        if dt_psi is not None:
            velocity_field = self.calc_vortex_velocity_field(dt_psi, self.psi[0])

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



    
