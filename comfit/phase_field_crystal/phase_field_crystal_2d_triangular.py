from comfit.phase_field_crystal.phase_field_crystal import PhaseFieldCrystal

import numpy as np
from tqdm import tqdm
from scipy.optimize import fsolve
import scipy as sp
import matplotlib.pyplot as plt
from pprint import pprint

from comfit.tool.tool_print_in_color import tool_print_in_color

class PhaseFieldCrystal2DTriangular(PhaseFieldCrystal):
    def __init__(self, nx, ny, **kwargs):
        """
        Initializes a phase field crystal system in 2D with a triangular crystal structure.

        Parameters
        ----------
        nx : int
            The number of unit cells in the x direction.
        ny : int
            The number of unit cells in the y direction.
        kwargs : dict
            Additional arguments to set as attributes, possibly overwriting default values.

        Returns
        -------
        None
        """

        # Type of the system
        self.type = 'PhaseFieldCrystal2DTriangular'
        self.dim = 2
        self.nx = nx
        self.ny = ny

        # Default simulation parameters
        self.micro_resolution = [7, 12]
        self.psi0 = -0.3
        self.r = -0.3
        self.t = 0
        self.v = 1
        self.dt = 0.1

        # If there are additional arguments provided, set them as attributes
        # Possibly overwriting the default values
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.xRes = nx * self.micro_resolution[0]
        self.yRes = ny * self.micro_resolution[1]

        a0 = 2 * np.pi * 2 / np.sqrt(3)
        self.a = a0 * np.array([[1, 0], [1 / 2, np.sqrt(3) / 2], [1 / 2, -np.sqrt(3) / 2]])
        self.q = np.array([[np.sqrt(3) / 2, -1 / 2], [0, 1], [-np.sqrt(3) / 2, -1 / 2]])

        # Set the number of reciprocal modes
        self.number_of_reciprocal_lattice_modes = 3
        self.number_of_primary_reciprocal_lattice_modes = 3

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = np.sqrt(3) * a0 / self.micro_resolution[1]

        self.type_of_evolution = kwargs.get('type_of_evolution', 'conserved')
        if self.type_of_evolution == 'unconserved':
            self.psi0_proto, self.A_proto = self.calc_proto_amplitudes_unconserved()
        else:
            if self.type_of_evolution != 'conserved':
                tool_print_in_color('Warning: type_of_evolution should be either conserved or unconserved. Setting to conserved.', 'red')
            self.A_proto = self.calc_proto_amplitudes_conserved()   

        # Set the elastic constants
        self.el_lambda_proto = 3 * self.A_proto ** 2
        self.el_mu_proto = 3 * self.A_proto ** 2
        self.el_gamma_proto = 0

        bool_is_for_properties_calculation = kwargs.get('for_properties_calculation', False)

        if not bool_is_for_properties_calculation:
            tool_print_in_color('Initiating a 2D triangular PFC model.', 'green')
            kwargs.pop('for_properties_calculation', None)
            pfc = PhaseFieldCrystal2DTriangular(1,1,for_properties_calculation=True, type_of_evolution=self.type_of_evolution, **kwargs)
            final_strain, self.psi0, self.A, self.el_lambda, self.el_mu, self.el_gamma = pfc.calc_strained_amplitudes()     

            # Calculate S1111_ref
            dxxxpsi = pfc.ifft(pfc.dif[0] * pfc.dif[0] * pfc.dif[0] * pfc.fft(pfc.psi)).real
            dxyypsi = pfc.ifft(pfc.dif[0] * pfc.dif[1] * pfc.dif[1] * pfc.fft(pfc.psi)).real
            dyyypsi = pfc.ifft(pfc.dif[1] * pfc.dif[1] * pfc.dif[1] * pfc.fft(pfc.psi)).real

            self.S111111_ref = np.mean(dxxxpsi ** 2)
            self.Skkllmm_ref = np.mean(dxxxpsi ** 2 + 3 * dxxxpsi*dxyypsi + 3 * dxyypsi**2 + dyyypsi**2)

        else:
            self.A = self.A_proto
            self.el_lambda = self.el_lambda_proto
            self.el_mu = self.el_mu_proto
            self.el_gamma = self.el_gamma_proto
            
        self.eta0 = np.array([self.A, self.A, self.A])

        # Initialize the BaseSystem
        kwargs['xRes'] = self.xRes
        kwargs['yRes'] = self.yRes
        kwargs['dx'] = self.dx
        kwargs['dy'] = self.dy
        kwargs['dt'] = self.dt
        super().__init__(self.dim,  **kwargs)
        
        # Set the a0
        self.a0 = a0
        self.bool_has_defined_length_scale = True
        
        if not bool_is_for_properties_calculation:
            self.conf_apply_distortion([[final_strain,0],[0,final_strain]], update_q_and_a_vectors=True)
            self.a0 = self.a0 *(1+final_strain)



    def calc_proto_amplitudes_conserved(self):
        """Calculate the proto-amplitudes for the system.

        Returns
        -------
        float
            The proto-amplitudes for the system.
        """

        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v
        
        # Testing which of the three amplitudes give the lowest free energy
        A = 0
        A_tmp = 1 / (15 * v) * (-t - 3 * v * psi0 - np.sqrt(t**2 - 15 * v * r - 24 * t * v * psi0 - 36 * v ** 2 * psi0 ** 2)) 
        if self.calc_free_energy_from_proto_amplitudes(psi0, A_tmp) < self.calc_free_energy_from_proto_amplitudes(psi0, A):
            A = A_tmp
        A_tmp = 1 / (15 * v) * (-t - 3 * v * psi0 + np.sqrt(t**2 - 15 * v * r - 24 * t * v * psi0 - 36 * v ** 2 * psi0 ** 2)) 
        if self.calc_free_energy_from_proto_amplitudes(psi0, A_tmp) < self.calc_free_energy_from_proto_amplitudes(psi0, A):
            A = A_tmp
        return A

    def calc_proto_amplitudes_unconserved(self):
        """Calculate the proto-amplitudes for the system.

        Returns
        -------
        tuple
            The proto-amplitudes for the system.
        """

        psi0 = 0
        A = 0
        free_energy = self.calc_free_energy_from_proto_amplitudes(psi0, A)

        for psi00 in np.linspace(-1, 1, 10):
            for A0 in np.linspace(0.01, 1, 5):
                [psi0_tmp, A_tmp] = fsolve(self.calc_proto_amplitude_equations_unconserved, [psi00, A0])
                free_energy_tmp = self.calc_free_energy_from_proto_amplitudes(psi0_tmp, A_tmp)

                if free_energy_tmp <= free_energy:
                    psi0 = psi0_tmp
                    A = A_tmp
                    free_energy = free_energy_tmp

        return psi0, A

    def calc_proto_amplitude_equations_unconserved(self, vars):
        """
        Calculates the equations for the proto-amplitudes for the system in the case of conserved dynamics.

        Parameters
        ----------
        vars : tuple
            The proto-amplitudes for the system.

        Returns
        -------
        list
            The equations for the proto-amplitudes for the system.
        """
        r = self.r
        t = self.t
        v = self.v
        
        psi0, A = vars
        eq1 = 12 * A**3 * self.v + 6 * A**2 * (self.t + 3 * self.v * psi0) + psi0 * (1 + self.r + self.t * psi0 + self.v * psi0**2)
        eq2 = self.r + 15 * A**2 * self.v + 2 * A * (self.t + 3 * self.v * psi0) + psi0 * (2 * self.t + 3 * self.v * psi0) # *A not necessary  

        return [eq1, eq2]

    def calc_free_energy_from_proto_amplitudes(self, psi0, A):
        """
        Calculates the free energy of the system from the proto-amplitudes.

        Parameters
        ----------
        psi0 : float
            The average value of psi.
        A : float
            The proto-amplitude.

        Returns
        -------
        float
            The free energy of the system.
        """
        r = self.r
        t = self.t
        v = self.v

        return np.pi**2 / (3 * np.sqrt(3)) * (270 * A**4 * v + 48 * A**3 * (t + 3 * v * psi0) + psi0**2 * (6 + 6 * r + 4 * t * psi0 + 3 * v * psi0**2) + 36 * A**2 * (r + psi0 * (2 * t + 3 * v * psi0)))

    def calc_L_f(self):
        """
        Calculates the L operator in Fourier space.

        Returns
        -------
        np.ndarray
            The L operator in Fourier space.
        """
        k2 = self.calc_k2()
        return 1 - k2

    def calc_L_sum_f(self):
        """
        Calculates the sum of the L operators in Fourier space. Needed for stress calculation functions.

        Returns
        -------
        int
            The L operator in Fourier space.
        """
        return 1