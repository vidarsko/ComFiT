from comfit.models.phase_field_crystal import PhaseFieldCrystal

import numpy as np
from tqdm import tqdm
from scipy.optimize import fsolve
import scipy as sp
import matplotlib.pyplot as plt
from pprint import pprint

class PhaseFieldCrystal1DPeriodic(PhaseFieldCrystal):
    def __init__(self, nx, **kwargs):
        """Initializes a phase field crystal system in 1D with a periodic crystal structure.

        Args:
            nx: The number of unit cells in the x direction.

        Returns:
            The system object representing the PhaseFieldCrystal1DPeriodic simulation.
        """

        # Type of the system
        self.type = 'PhaseFieldCrystal1DPeriodic'
        self.dim = 1
        self.nx = nx

        # Default simulation parameters
        self.micro_resolution = kwargs.get('micro_resolution',[5])
        self.psi0 = -0.3
        self.r = -0.3
        self.t = 0
        self.v = 1
        self.dt = 0.1

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.xRes = nx * self.micro_resolution[0]

        a0 = 2 * np.pi
        self.a = a0 * np.array([[1]])
        self.q = np.array([[1]])

        # Set the number of reciprocal modes
        self.number_of_reciprocal_lattice_modes = 1
        self.number_of_primary_reciprocal_lattice_modes = 1

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]

        self.type_of_evolution = kwargs.get('type_of_evolution', 'conserved')
        self.A = self.calc_proto_amplitudes_conserved()
        self.eta0 = np.array([self.A])

        # Set the elastic constants
        self.el_mu = 3 * self.A ** 2
        self.el_lambda = 3 * self.A ** 2
        self.el_gamma = 0

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, 
                         dx=self.dx,  dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

    def calc_proto_amplitudes_conserved(self):
        """Calculates the proto-amplitudes for the system. 

        Args:
            None
        
        Returns:
            The proto-amplitudes for the system.
        """


        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v

        A = 0
        
        A_tmp = 1/np.sqrt(3*v)*np.sqrt(- r - 2*t*psi0 - 3*v*psi0**2)
        if self.calc_free_energy_from_proto_amplitudes(psi0, A_tmp) < self.calc_free_energy_from_proto_amplitudes(psi0, A):
            A = A_tmp

        return A

    def calc_free_energy_from_proto_amplitudes(self, psi0, A):
        """Calculates the free energy of the system from the proto-amplitudes.

        Args:
            psi0: The average value of psi.
            A: The proto-amplitude.

        Returns:
            The free energy of the system.
        """

        r = self.r
        t = self.t
        v = self.v

        return 2*np.pi*(1/2*psi0**2 + 1/2*r*(psi0**2 + 2*A**2) + 1/3*t*(psi0**3 + 6*psi0*A**2) + 1/4*v*(psi0**4 + 12*psi0**2*A**2 + 6*A**4))

    def calc_L_f(self):
        """Calculates the L operator in Fourier space.

        Args:
            None
        Returns:
            The L operator in Fourier space.
        """
        k2 = self.calc_k2()
        return 1 - k2

    def calc_omega_f(self):
        """Calculates the free energy of the system.

        Args:
            None
            
        Returns:
            The free energy of the system.
        """
        
        k2 = self.calc_k2()
        return - k2 * (self.r + (1 - k2)**2)