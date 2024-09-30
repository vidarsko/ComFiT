from comfit.models.phase_field_crystal import PhaseFieldCrystal

import numpy as np
from tqdm import tqdm
from scipy.optimize import fsolve
import scipy as sp
import matplotlib.pyplot as plt
from pprint import pprint

class PhaseFieldCrystal3DBodyCenteredCubic(PhaseFieldCrystal):
    def __init__(self, nx, ny, nz, **kwargs):
        """Initializes a phase field crystal system in 3D with a body centered cubic crystal structure.

        Args:
            nx: The number of unit cells in the x direction.
            ny: The number of unit cells in the y direction.
            nz: The number of unit cells in the z direction.

        Returns:
            The system object representing the PhaseFieldCrystal3DBodyCenteredCubic simulation.
        """

        # Type of the system
        self.type = 'PhaseFieldCrystal3DBodyCenteredCubic'
        self.dim = 3
        self.nx = nx
        self.ny = ny
        self.nz = nz

        # Default simulation parameters
        self.micro_resolution = [7, 7, 7]
        self.psi0 = -0.325
        self.r = -0.3
        self.t = 0
        self.v = 1
        self.dt = 0.1

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.xRes = nx * self.micro_resolution[0]
        self.yRes = ny * self.micro_resolution[1]
        self.zRes = nz * self.micro_resolution[2]

        a0 = 2 * np.pi * np.sqrt(2)
        self.a = a0 / 2 * np.array([[-1, 1, 1],
                                    [1, -1, 1],
                                    [1, 1, -1],
                                    [1, 1, 1]])

        self.q = np.array([[0, 1, 1],
                           [1, 0, 1],
                           [1, 1, 0],
                           [0, -1, 1],
                           [-1, 0, 1],
                           [-1, 1, 0]]) / np.sqrt(2)

        # Set the number of reciprocal modes
        self.number_of_reciprocal_lattice_modes = 6
        self.number_of_primary_reciprocal_lattice_modes = 6

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]
        self.dz = a0 / self.micro_resolution[2]

        self.type_of_evolution = kwargs.get('type_of_evolution', 'conserved')
        self.A = self.calc_proto_amplitudes_conserved()
        self.eta0 = np.array([self.A, self.A, self.A, self.A, self.A, self.A])

        # Set the elastic constants
        self.el_lambda = 4 * self.A ** 2
        self.el_mu = 4 * self.A ** 2
        self.el_gamma = - 4*self.A**2

        if not bool_is_for_properties_calculation:
            tool_print_in_color('Initiating a 2D triangular PFC model.', 'green')
            pfc = PhaseFieldCrystal2DTriangular(1,1,for_properties_calculation=True, type_of_evolution=self.type_of_evolution)
            final_strain, self.psi0, self.A, self.el_lambda, self.el_mu, self.el_gamma = pfc.calc_strained_amplitudes()  

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes, zRes=self.zRes,
                         dx=self.dx, dy=self.dy, dz=self.dz, dt=self.dt)

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

        # Testing which of the three amplitudes give the lowest free energy
        A=0
        A_tmp = 1/(45*v)*( -2*t - 6*v*psi0 + np.sqrt(4*t**2 - 45*r*v - 66*t*v*psi0 - 99*v**2*psi0**2))
        if self.calc_free_energy_from_proto_amplitudes(psi0, A_tmp) < self.calc_free_energy_from_proto_amplitudes(psi0, A):
            A = A_tmp
        A_tmp = 1/(45*v)*( -2*t - 6*v*psi0 - np.sqrt(4*t**2 - 45*r*v - 66*t*v*psi0 - 99*v**2*psi0**2))
        if self.calc_free_energy_from_proto_amplitudes(psi0, A_tmp) < self.calc_free_energy_from_proto_amplitudes(psi0, A):
            A = A_tmp
        return A

    def calc_free_energy_from_proto_amplitudes(self, psi0, A):
        """Calculates the free energy of the phase-field crystal from the proto amplitudes.

        Args:
            psi0: The average value of psi.
            A: The proto-amplitude.

        Returns:
            The free energy of the phase-field crystal.
        """
        r = self.r
        t = self.t
        v = self.v

        return 4*np.sqrt(2)*np.pi**3/3*(1620*A**4*v + 192*A**3*(t + 3*v*psi0) + psi0**2*(6 + 6*r + 4*t*psi0 + 3*v*psi0**2) + 72*A**2*(r + psi0*(2*t + 3*v*psi0)))

    def calc_L_f(self):
        """Calculates the L operator in Fourier space.

        Args:
            None
        Returns:
            The L operator in Fourier space.
        """
        k2 = self.calc_k2()
        return 1 - k2

    def calc_L_sum_f(self):
        """Calculates the sum of the L operators in Fourier space. Needed for stress calculation functions.

        Args:
            None
        Returns:
            The L operator in Fourier space.
        """
        return 1