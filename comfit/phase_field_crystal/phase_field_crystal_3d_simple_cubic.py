from comfit.phase_field_crystal.phase_field_crystal import PhaseFieldCrystal

import numpy as np
from tqdm import tqdm
from scipy.optimize import fsolve
import scipy as sp
import matplotlib.pyplot as plt
from pprint import pprint

from comfit.tool.tool_print_in_color import tool_print_in_color

class PhaseFieldCrystal3DSimpleCubic(PhaseFieldCrystal):
    def __init__(self, nx, ny, nz, **kwargs):
        """
        Initializes a phase field crystal system in 3D with a simple cubic crystal structure.

        Parameters
        ----------
        nx : int
            The number of unit cells in the x direction.
        ny : int
            The number of unit cells in the y direction.
        nz : int
            The number of unit cells in the z direction.
        \*\*kwargs : dict
            Additional arguments to set as attributes.

        Returns
        -------
        None
        """
        # Type of the system
        self.type = 'PhaseFieldCrystal3DSimpleCubic'
        self.dim = 3
        self.nx = nx
        self.ny = ny
        self.nz = nz

        # Default simulation parameters
        self.micro_resolution = [5, 5, 5]
        self.psi0 = -0.325
        self.r = -0.3
        self.t = 0
        self.v = 1
        self.dt = 0.1
        # TODO: Some of this code seems to be reprinted. Maybe it should be moved to the BaseSystem class? (Vidar 18.12.23)

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.xRes = nx * self.micro_resolution[0]
        self.yRes = ny * self.micro_resolution[1]
        self.zRes = nz * self.micro_resolution[2]

        a0 = 2 * np.pi
        self.a = a0 * np.array([[1, 0, 0],
                                [0, 1, 0],
                                [0, 0, 1]])

        self.q = np.array([[1, 0, 0],
                           [0, 1, 0],
                           [0, 0, 1],
                           [0, 1, 1],
                           [1, 0, 1],
                           [1, 1, 0],
                           [0, -1, 1],
                           [-1, 0, 1],
                           [-1, 1, 0],
                           [-1, 1, 1],
                           [1, -1, 1],
                           [1, 1, -1],
                           [1, 1, 1]], dtype=float)

        # Set the number of reciprocal modes
        self.number_of_reciprocal_lattice_modes = 13
        self.number_of_primary_reciprocal_lattice_modes = 3

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]
        self.dz = a0 / self.micro_resolution[2]

        self.type_of_evolution = kwargs.get('type_of_evolution', 'conserved')
        self.A_proto, self.B_proto, self.C_proto = self.calc_proto_amplitudes_conserved()

        # Set the elastic constants
        self.el_lambda_proto = 16 * self.B_proto ** 2 + 128 * self.C_proto ** 2
        self.el_mu_proto = 16 * self.B_proto ** 2 + 128 * self.C_proto ** 2
        self.el_gamma_proto = 32*self.A_proto**2 - 16*self.B_proto**2 - 256*self.C_proto**2

        bool_is_for_properties_calculation = kwargs.get('for_properties_calculation', False)

        if not bool_is_for_properties_calculation:
            tool_print_in_color('Initiating a 3D simple cubic PFC model.', 'green')
            kwargs.pop('for_properties_calculation', None)
            pfc = PhaseFieldCrystal3DSimpleCubic(1,1,1,for_properties_calculation=True, type_of_evolution=self.type_of_evolution, **kwargs)
            final_strain, self.psi0, self.A, self.B, self.C, self.el_lambda, self.el_mu, self.el_gamma = pfc.calc_strained_amplitudes()     
        else:
            self.A = self.A_proto
            self.B = self.B_proto
            self.C = self.C_proto
            self.el_lambda = self.el_lambda_proto
            self.el_mu = self.el_mu_proto
            self.el_gamma = self.el_gamma_proto


        self.eta0 = np.array([self.A, self.A, self.A,
                     self.B, self.B, self.B, self.B, self.B, self.B,
                     self.C, self.C, self.C, self.C])

        # Initialize the BaseSystem
        kwargs['xRes'] = self.xRes
        kwargs['yRes'] = self.yRes
        kwargs['zRes'] = self.zRes
        kwargs['dx'] = self.dx
        kwargs['dy'] = self.dy
        kwargs['dz'] = self.dz
        kwargs['dt'] = self.dt
        super().__init__(self.dim,  **kwargs)

        # Set the a0
        self.a0 = a0
        self.bool_has_defined_length_scale = True

        if not bool_is_for_properties_calculation:
            self.conf_apply_distortion([[final_strain,0,0],[0,final_strain,0],[0,0,final_strain]], update_q_and_a_vectors=True)
            self.a0 = self.a0 * (1+final_strain)

    def calc_free_energy_from_proto_amplitudes(self, psi0, A, B, C):
        """
        Calculates the free energy of the system from the proto-amplitudes.

        Parameters
        ----------
        psi0 : float
            The value of psi0.
        A : float
            The value of A.
        B : float
            The value of B.
        C : float
            The value of C.

        Returns
        -------
        float
            The free energy of the system.
        """
        r = self.r
        t = self.t
        v = self.v
        return 2*np.pi**3/3*(48*C**2*r + 270*A**4*v + 1620*B**4*v + 576*A**3*C*v + 648*C**4*v + 96*C**2*t*psi0 + 6*(36 + r + 24*C**2*v)*psi0**2 + 4*t*psi0**3 + 3*v*psi0**4 + 192*B**3*(t + 3*v*psi0) + 576*A*B*C*(t + 3*v*(3*B + psi0)) + 36*A**2*( r + 96*B**2*v + 36*C**2*v + 2*t*psi0 + 3*v*psi0**2 + 8*B*( t + 3*v*psi0)) + 72*B**2*(r + 54*C**2*v + psi0*(2*t + 3*v*psi0)) )

    def calc_proto_amplitudes_conserved(self):
        """
        Calculates the proto-amplitudes for the system.

        Parameters
        ----------
        None

        Returns
        -------
        tuple
            The proto-amplitudes for the system.
        """
        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v

        A = 0
        B = 0
        C = 0
        free_energy = self.calc_free_energy_from_proto_amplitudes(psi0, A, B, C)

        for A0 in np.linspace(0.01, 1, 10):
            for B0 in [A0/2]:
                for C0 in [A0/4]:
                    [A_tmp, B_tmp, C_tmp] = fsolve(self.calc_proto_amplitude_equations_conserved, [A0, B0, C0])
                    free_energy_tmp = self.calc_free_energy_from_proto_amplitudes(psi0, A_tmp, B_tmp, C_tmp)

                    if free_energy_tmp <= free_energy:
                        A = A_tmp
                        B = B_tmp
                        C = C_tmp
                        free_energy = free_energy_tmp

        return A, B, C

    def calc_proto_amplitude_equations_conserved(self, vars):
        """
        Calculates the equations for the proto-amplitudes for the system in the case of conserved dynamics.

        Parameters
        ----------
        vars : list
            The proto-amplitudes for the system.

        Returns
        -------
        list
            The equations for the proto-amplitudes for the system in the case of conserved dynamics.
        """
        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v

        A, B, C = vars

        eq1 = 15*A**3*v + 24*A**2*C*v + 8*B*C*(t + 3*v*(3*B + psi0)) + A*( r + 96*B**2*v + 36*C**2*v + 2*t*psi0 + 3*v*psi0**2 + 8*B*(t + 3*v*psi0))
        eq2 = 45*B**3*v + 4*B**2*(t + 3*v*psi0) + 2*A*(A + 2*C)*(t + 3*v*psi0) + B*( r + 48*A**2*v + 72*A*C*v + 54*C**2*v + 2*t*psi0 + 3*v*psi0**2)
        eq3 = 27*C**3*v + C*( r + 27*A**2*v + 81*B**2*v + 2*t*psi0 + 3*v*psi0**2) + 6*A*(A**2*v + 9*B**2*v + B*(t + 3*v*psi0))

        return [eq1, eq2, eq3]

    def calc_L_f(self):
        """
        Calculates the L operator in Fourier space.

        Parameters
        ----------
        None

        Returns
        -------
        numpy.ndarray
            The L operator in Fourier space.
        """
        k2 = self.calc_k2()
        return (1 - k2)*(2 - k2)*(3 - k2)

    def calc_L_sum_f(self):
        """
        Calculates the sum of the L operators in Fourier space. Needed for stress calculation functions.

        Parameters
        ----------
        None

        Returns
        -------
        numpy.ndarray
            The L operator in Fourier space.
        """
        k2 = self.calc_k2()
        return 11 - 12*k2 + 3*k2**2