from typing import Optional, Union, Literal, TYPE_CHECKING, List, Dict, Any

# General packages
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

# Local packages
from comfit.core import BaseSystem
from comfit.tool import tool_colormap

class QuantumMechanics(BaseSystem):
    def __init__(self : 'QuantumMechanics', 
            dimension: int, 
            **kwargs: Dict[str, Any]
            ) -> None:
        """Initializes a quamtum mechanics system evolving according to the Schrödinger equation

        Parameters
        ----------
        dimension : int
            The dimension of the system.
        kwargs : dict, optional
            Optional keyword arguments to set additional parameters. 
            See https://comfitlib.com/ClassBaseSystem/.

        Examples
        --------
        >>> qm = cf.QuantumMechanics(3,xRes=101,yRes=101,zRes=101)
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

    def __str__(
        self: 'QuantumMechanics',
        ) -> str:
        """Returns a string representation of the QuantumMechanics system.

        This method is called when str() is invoked on an instance of the class.

        Returns
        -------
        str: String representation 'QuantumMechanics' of the system
        """
        return 'QuantumMechanics'

    def conf_initial_condition_Gaussian(
            self : 'QuantumMechanics',
            position: Optional[Union[List[float], np.ndarray]] = None,
            width: Optional[float] = None,
            initial_velocity: Optional[Union[float, list[float]]] = None
            ) -> None:
        """Set the initial condition to a Gaussian wavepacket.

        Parameters
        ----------
        position : array_like, optional
            The position of the wavepacket.
        width : float, optional
            The width of the wavepacket.
        initial_velocity : float or array_like, optional
            The initial velocity of the wavepacket.

        Returns
        -------
        None
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
                self.psi = self.psi * np.exp(1j * (v0[0] * self.x + v0[1] * self.y))

            elif self.dim == 3:
                v0 = initial_velocity
                self.psi = self.psi * np.exp(1j * (v0[0] * self.x + v0[1] * self.y + v0[2] * self.z))
        
        self.psi_f = sp.fft.fftn(self.psi)


    def conf_harmonic_potential(
            self : 'QuantumMechanics', 
            trapping_strength: Optional[float] = None
            ) -> None:
        """Set the external potential to a harmonic trap with R_tf being the thomas fermi radius
        
        Parameters
        ----------
        trapping_strength : float, optional
            The strength of the trapping potential.

        Returns
        -------
        None       
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

    def conf_hydrogen_state(
            self: 'QuantumMechanics', 
            n: int, 
            l: int, 
            m: int
            ) -> None:
        """Set the initial condition to a hydrogen state with quantum numbers n,l,m

        Parameters
        ----------
        n : int 
            principal quantum number
        l : int
            azimuthal quantum number
        m : int 
            magnetic quantum number

        Returns
        -------
        None
            Configures self.psi and self.psi_f.
        """

        self.psi = self.calc_hydrogen_state(n,l,m)
        self.psi_f = sp.fft.fftn(self.psi)

    def conf_wavefunction(
            self: 'QuantumMechanics', 
            psi: np.ndarray
            ) -> None:
        """Set the initial condition to a given wavefunction.

        Parameters
        ----------
        psi : np.ndarray
            The wavefunction.

        Returns
        -------
        None
            Configures self.psi and self.psi_f.
        """
        self.psi = psi
        self.psi_f = sp.fft.fftn(self.psi)   

    def calc_hydrogen_state(
            self: 'QuantumMechanics', 
            n: int, 
            l: int, 
            m: int
            ) -> np.ndarray:
        """Calculate the hydrogen state with quantum numbers n,l,m

        Parameters
        ----------
        n : int
            principal quantum number
        l : int
            azimuthal quantum number
        m : int
            magnetic quantum number

        Returns
        -------
        np.ndarray
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
    def calc_nonlinear_evolution_function_f(
            self: 'QuantumMechanics', 
            psi: np.ndarray, 
            t: float
            ) -> np.ndarray:
        """Calculate the nonlinear evolution function f of the Schrödinger equation

        Parameters
        ----------
        psi : np.ndarray 
            The wavefunction.
        t : float
            Time
        
        Returns
        -------
        np.ndarray
            The nonlinear evolution function f.
        """

        if callable(self.V_ext):
            potential = self.V_ext(t)
        else:
            potential = self.V_ext
        return sp.fft.fftn((1j) * (-potential) * psi)


    ## Time evolution functions
    def evolve_schrodinger(
            self: 'QuantumMechanics', 
            number_of_steps: int, 
            method: Literal['ETD2RK', 'ETD4RK'] = 'ETD4RK'
            ) -> None:
        """Evolve the system according to the Schrödinger equation

        Parameters
        ----------
        number_of_steps : int
            Number of time steps to evolve the system.
        method : str
            The method to use for the time evolution. 
        
        Returns
        -------
        None
            Updates self.psi and self.psi_f.
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