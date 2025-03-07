from typing import Optional, Union, Literal

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

from comfit.core.base_system import BaseSystem
from comfit.tool.tool_colormap import tool_colormap

class QuantumMechanics(BaseSystem):
    def __init__(self, dimension: int, **kwargs) -> None:
        """Initializes a quamtum mechanics system evolving according to the Schrödinger equation

        Args:
            dimension: The dimension of the system.
            **kwargs : dict, optional. Optional keyword arguments to set additional parameters. Same as BaseSystem
                
        Returns:
            The system object representing the QuantumMechanics simulation.

        Example:
            qm = cf.QuantumMechanics(3,xRes=101,yRes=101,zRes=101)
            Creates a BoseEinsteinCondensate system with 3 dimensions and a spatial resolution of 101.
        """

        # First initialize the BaseSystem class
        super().__init__(dimension, **kwargs)

        # Type of the system
        self.psi = None
        self.psi_f = None
        self.type = 'QuantumMechanics1' #Quantum mechanical system with one particle.

        # Default simulation parameters
        self.V_ext = 0  # External potential

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __str__(self) -> str:
        """Output a string representation of the system

        Args:
            None

        Returns:
            A string representation of the system
        """
        return 'QuantumMechanics'

    def conf_initial_condition_Gaussian(
            self,
            position: Optional[list[complex]] = None,
            width: Optional[float] = None,
            initial_velocity: Optional[Union[float, list[float]]] = None
        ):
        """Set the initial condition to a Gaussian wavepacket.

        Args:
            position: The position of the wavepacket.
            width: The width of the wavepacket.
            initial_velocity: The initial velocity of the wavepacket.

        Returns:
            Configures self.psi and self.psi_f.
        """

        if position is None:
            position = self.rmid
        if width is None:
            width = self.size_x/10

        self.psi = np.sqrt(self.calc_Gaussian(position=position,width=width))
        
        if initial_velocity is not None:
            if self.dim == 1:
                v0 = initial_velocity
                self.psi = self.psi * np.exp(1j * v0 * self.x)

            elif self.dim == 2:
                v0 = initial_velocity
                self.psi = self.psi * np.exp(1j * (v0[0] * self.x + v0[1]*self.y))

            elif self.dim == 3:
                v0 = initial_velocity
                self.psi = self.psi * np.exp(1j * (v0[0] * self.x + v0[1]*self.y + v0[2]*self.z))
        
        self.psi_f = sp.fft.fftn(self.psi)


    def conf_harmonic_potential(self, trapping_strength: Optional[float] = None) -> None:
        """Set the external potential to a harmonic trap with R_tf being the thomas fermi radius
        
        Args:
            trapping_strength: The strength of the trapping potential.

        Returns:
            Configures self.V_ext.
        """

        if trapping_strength is None:
            trapping_strength = 1/(0.25*self.size_x)**2

        if self.dim == 1:
            return  trapping_strength*(self.x - self.xmid )**2

        if self.dim == 2:
            return trapping_strength*(((self.x-self.xmid)**2) +((self.y-self.ymid)**2) )
        if self.dim == 3:
            return trapping_strength * (((self.x - self.xmid) ** 2)
                                           + ((self.y - self.ymid) ** 2)
                                           +((self.z - self.zmid) ** 2) )

    def conf_hydrogen_state(self, n: int, l: int, m: int) -> None:
        """Set the initial condition to a hydrogen state with quantum numbers n,l,m

        Args:
            n: principal quantum number
            l: azimuthal quantum number
            m: magnetic quantum number

        Returns:
            Configures self.psi and self.psi_f.
        """

        self.psi = self.calc_hydrogen_state(n,l,m)
        self.psi_f = sp.fft.fftn(self.psi)

    def conf_wavefunction(self, psi: np.ndarray) -> None:
        """Set the initial condition to a given wavefunction.

        Args:
            psi: The wavefunction.

        Returns:
            Configures self.psi and self.psi_f.
        """
        self.psi = psi
        self.psi_f = sp.fft.fftn(self.psi)   

    def calc_hydrogen_state(self, n: int, l: int, m: int) -> np.ndarray:
        """Calculates the hydrogen state with quantum numbers n,l,m

        Args:
            n: principal quantum number
            l: azimuthal quantum number
            m: magnetic quantum number

        Returns:
            Wavefunction of the hydrogen state. (np.ndarray)
        """

        if not self.dim == 3:
            raise(Exception('The hydrogen state is only defined in 3D'))

        r = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        theta = np.arccos(self.z/r)
        theta[np.isnan(theta)] = 0
        phi = np.arctan2(self.y,self.x)

        R = (2*r/n)**l * np.exp(-r/n) *  sp.special.genlaguerre(n-l-1,2*l+1)(2*r/n)
        Y = sp.special.sph_harm(m,l,phi,theta)
        
        return R*Y
    

    ## Calculation functions
    def calc_nonlinear_evolution_function_f(self, psi: np.ndarray, t: float) -> np.ndarray:
        """Calculates the nonlinear evolution function f of the Schrödinger equation

        Args:
            psi: The wavefunction.
            t: Time.
        
        Returns:
            The nonlinear evolution function f.
        """

        if callable(self.V_ext):
            potential = self.V_ext(t)
        else:
            potential = self.V_ext
        return sp.fft.fftn((1j) * (-potential) * psi)


    ## Time evolution functions
    def evolve_schrodinger(self, number_of_steps: int, method: Literal['ETD2RK', 'ETD4RK'] = 'ETD4RK') -> None:
        """Evolve the system according to the Schrödinger equation

        Args:
            number_of_steps: Number of time steps to evolve the system.
            method: The method to use for the time evolution. 
        
        Returns:
            None: Does not return anything. Just evolves the system according to the Schrödinger equation.
        """

        omega_f = -1j / 2 * self.calc_k2()
        if method == 'ETD2RK':
            integrating_factors_f = self.calc_evolution_integrating_factors_ETD2RK(omega_f)
            solver = self.evolve_ETD2RK_loop
        elif method == 'ETD4RK':
            integrating_factors_f = self.calc_evolution_integrating_factors_ETD4RK(omega_f)
            solver = self.evolve_ETD4RK_loop
        else:
            raise(Exception('This method is not implemented'))
        for n in range(number_of_steps):
            self.psi, self.psi_f = solver(integrating_factors_f, self.calc_nonlinear_evolution_function_f,
                                                           self.psi, self.psi_f)