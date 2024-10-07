from comfit.models.phase_field_crystal import PhaseFieldCrystal

import numpy as np
from tqdm import tqdm
from scipy.optimize import fsolve
import scipy as sp
import matplotlib.pyplot as plt
from pprint import pprint

from comfit.tool.tool_print_in_color import tool_print_in_color

class PhaseFieldCrystal3DFaceCenteredCubic(PhaseFieldCrystal):
    def __init__(self, nx, ny, nz, **kwargs):
        """Initializes a phase field crystal system in 3D with a face centered cubic crystal structure.

        Args:
            nx: The number of unit cells in the x direction.
            ny: The number of unit cells in the y direction.
            nz: The number of unit cells in the z direction.

        Returns:
            The system object representing the PhaseFieldCrystal3DFaceCenteredCubic simulation.
        """

        # Type of the system
        self.type = 'PhaseFieldCrystal3DFaceCenteredCubic'
        self.dim = 3
        self.nx = nx
        self.ny = ny
        self.nz = nz

        # Default simulation parameters
        self.micro_resolution = [11, 11, 11]
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

        a0 = 2 * np.pi * np.sqrt(3)
        self.a = a0 / 2 * np.array([[0, 1, 1],
                                    [1, 0, 1],
                                    [1, 1, 0],
                                    [0, -1, 1],
                                    [-1, 0, 1],
                                    [-1, 1, 0]], dtype=float)

        self.q = np.array([[-1, 1, 1],
                           [1, -1, 1],
                           [1, 1, -1],
                           [1, 1, 1],
                           [2, 0, 0],
                           [0, 2, 0],
                           [0, 0, 2]], dtype=float) / np.sqrt(3)

        # Set the number of reciprocal modes
        self.number_of_reciprocal_lattice_modes = 7
        self.number_of_primary_reciprocal_lattice_modes = 4

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]
        self.dz = a0 / self.micro_resolution[2]

        self.type_of_evolution = kwargs.get('type_of_evolution', 'conserved')
        self.A_proto, self.B_proto = self.calc_proto_amplitudes_conserved()
    
        # Set the elastic constants
        self.el_lambda_proto = 32/81 * self.A_proto ** 2
        self.el_mu_proto = 32/81 * self.A_proto ** 2
        self.el_gamma_proto = 32/81 * (2*self.B_proto**2 - self.A_proto**2)

        bool_is_for_properties_calculation = kwargs.get('for_properties_calculation', False)
        if not bool_is_for_properties_calculation:
            tool_print_in_color('Initiating a 3D fcc PFC model.', 'green')
            kwargs.pop('for_properties_calculation', None)
            pfc = PhaseFieldCrystal3DFaceCenteredCubic(1,1,1,for_properties_calculation=True, type_of_evolution=self.type_of_evolution, **kwargs)
            final_strain, self.psi0, self.A, self.B, self.el_lambda, self.el_mu, self.el_gamma = pfc.calc_strained_amplitudes()  
        else:
            self.A = self.A_proto
            self.B = self.B_proto
            self.el_lambda = self.el_lambda_proto
            self.el_mu = self.el_mu_proto
            self.el_gamma = self.el_gamma_proto
            
        self.eta0 = np.array([self.A, self.A, self.A, self.A, self.B, self.B, self.B])

        # Initialize the BaseSystem
        kwargs['xRes'] = self.xRes
        kwargs['yRes'] = self.yRes
        kwargs['zRes'] = self.zRes
        kwargs['dx'] = self.dx
        kwargs['dy'] = self.dy
        kwargs['dz'] = self.dz
        kwargs['dt'] = self.dt
        super().__init__(self.dim, **kwargs)

        # Set the a0
        self.a0 = a0
        self.bool_has_defined_length_scale = True

        if not bool_is_for_properties_calculation:
            self.conf_apply_distortion([[final_strain,0,0],[0,final_strain,0],[0,0,final_strain]], update_q_and_a_vectors=True)
            self.a0 = self.a0 * (1+final_strain)

        
    def calc_proto_amplitudes_conserved(self):
        """Calculates the proto-amplitudes for the system.

        Args:
            None
        
        Returns:
            The proto-amplitudes for the system.
        """

        A = 0
        B = 0
        free_energy = self.calc_free_energy_from_proto_amplitudes(self.psi0, A, B)

        for A0 in np.linspace(0.01, 1, 10):
            for B0 in [A0/2]:
                [A_tmp, B_tmp] = fsolve(self.calc_proto_amplitude_equations_conserved, [A0, B0])
                free_energy_tmp = self.calc_free_energy_from_proto_amplitudes(self.psi0, A_tmp, B_tmp)

                if free_energy_tmp <= free_energy:
                    A = A_tmp
                    B = B_tmp
                    free_energy = free_energy_tmp

        return A, B

    def calc_proto_amplitude_equations_conserved(self, vars):
        """Calculates the equations for the proto-amplitudes for the system in the case of conserved dynamics

        Args:
            vars: The proto-amplitudes for the system.
        
        Returns:
            The equations for the proto-amplitudes for the system in the case of conserved dynamics.
        """
        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v

        A, B = vars
        eq1 = A*( r + 27*A**2*v + 36*B**2*v + 2*t*psi0 + 3*v*psi0**2 + 6*B*(t + 3*v*psi0))
        eq2 = 15*B**3*v + 4*A**2*(t + 3*v*psi0) + B*( r + 48*A**2*v + 2*t*psi0 + 3*v*psi0**2)

        return [eq1, eq2]

    def calc_free_energy_from_proto_amplitudes(self, psi0, A, B):
        """Calculates the free energy of the system from the proto-amplitudes.

        Args:
            psi0: The average value of psi.
            A: The proto-amplitude.
            B: The proto-amplitude.

        Returns:
            The free energy of the system.
        """

        r = self.r
        t = self.t
        v = self.v
        return 2*np.pi**3/np.sqrt(3)*(1944*A**4*v + 810*B**4*v + psi0**2*(32 + 18*r + 12*t*psi0 + 9*v*psi0**2) + 108*B**2*(r + psi0*(2*t + 3*v*psi0)) + 144*A**2*( r + 36*B**2*v + 6*B*(t + 3*v*psi0) + psi0*(2*t + 3*v*psi0)))

    def calc_L_f(self):
        """Calculates the L operator in Fourier space.

        Args:
            None
        Returns:
            The L operator in Fourier space.
        """
        k2 = self.calc_k2()
        return (1 - k2)*(4/3-k2)

    def calc_L_sum_f(self):
        """Calculates the sum of the L operators in Fourier space. Needed for stress calculation functions.

        Args:
            None
        Returns:
            The L operator in Fourier space.
        """
        return 7/3 - 2*k2