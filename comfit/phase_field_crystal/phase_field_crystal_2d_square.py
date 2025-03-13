from comfit.phase_field_crystal.phase_field_crystal import PhaseFieldCrystal

import numpy as np
from tqdm import tqdm
from scipy.optimize import fsolve
import scipy as sp
import matplotlib.pyplot as plt
from pprint import pprint

from comfit.tool.tool_print_in_color import tool_print_in_color

class PhaseFieldCrystal2DSquare(PhaseFieldCrystal):
    def __init__(self, nx, ny, **kwargs):
        """
        Initialize a phase field crystal system in 2D with a square crystal structure.

        Parameters
        ----------
        nx : int
            The number of unit cells in the x direction.
        ny : int
            The number of unit cells in the y direction.
        \*\*kwargs : dict, optional
            Additional keyword arguments to customize the simulation:
                micro_resolution : list
                    Resolution within each unit cell [x, y]
                psi0 : float
                    Average value of the density field
                r : float
                    Temperature parameter
                t : float
                    Parameter related to three-point correlation
                v : float
                    Parameter related to four-point correlation
                dt : float
                    Time step for simulation
                type_of_evolution : str
                    Type of dynamics ('conserved' or other)
                for_properties_calculation : bool
                    Whether this instance is for properties calculation

        Returns
        -------
        PhaseFieldCrystal2DSquare
            The system object representing the simulation.
        """

        # Type of the system
        self.type = 'PhaseFieldCrystal2DSquare'
        self.dim = 2
        self.nx = nx
        self.ny = ny

        # Default simulation parameters
        self.micro_resolution = [7, 7]
        self.psi0 = -0.3
        self.r = -0.3
        self.t = 0
        self.v = 1
        self.dt = 0.1

        # If there are additional arguments provided, set them as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.xRes = nx * self.micro_resolution[0]
        self.yRes = ny * self.micro_resolution[1]

        a0 = 2 * np.pi
        self.a = a0 * np.array([[1, 0], [0, 1]])
        self.q = np.array([[1, 0], [0, 1], [1, -1], [1, 1]], dtype=float)

        # Set the number of reciprocal modes
        self.number_of_reciprocal_lattice_modes = 4
        self.number_of_primary_reciprocal_lattice_modes = 2

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]

        self.type_of_evolution = kwargs.get('type_of_evolution', 'conserved')
        self.A_proto, self.B_proto = self.calc_proto_amplitudes_conserved()

        # Set the elastic constants
        self.el_lambda_proto = 16 * self.B_proto ** 2
        self.el_mu_proto = 16 * self.B_proto ** 2
        self.el_gamma_proto = 8 * self.A_proto ** 2 - 32 * self.B_proto ** 2

        bool_is_for_properties_calculation = kwargs.get('for_properties_calculation', False)

        if not bool_is_for_properties_calculation:
            tool_print_in_color('Initiating a 2D square PFC model.', 'green')
            kwargs.pop('for_properties_calculation', None)
            pfc = PhaseFieldCrystal2DSquare(1,1,for_properties_calculation=True, type_of_evolution=self.type_of_evolution, **kwargs)
            final_strain, self.psi0, self.A, self.B, self.el_lambda, self.el_mu, self.el_gamma = pfc.calc_strained_amplitudes()

            # Calculate S1111_ref
            dxxpsi = pfc.ifft(pfc.dif[0] * pfc.dif[0] * pfc.fft(pfc.psi)).real
            dxypsi = pfc.ifft(pfc.dif[0] * pfc.dif[1] * pfc.fft(pfc.psi)).real
            dyypsi = pfc.ifft(pfc.dif[1] * pfc.dif[1] * pfc.fft(pfc.psi)).real

            self.S1111_ref = np.mean(dxxpsi ** 2)
            self.Skkll_ref = np.mean(dxxpsi ** 2 + 2 * dxypsi ** 2 + dyypsi ** 2 )
            
        else:
            self.A = self.A_proto
            self.B = self.B_proto
            self.el_lambda = self.el_lambda_proto
            self.el_mu = self.el_mu_proto
            self.el_gamma = self.el_gamma_proto

        self.eta0 = np.array([self.A, self.A, self.B, self.B])

        # Initialize the BaseSystem
        kwargs['xRes'] = self.xRes
        kwargs['yRes'] = self.yRes
        kwargs['dx'] = self.dx
        kwargs['dy'] = self.dy
        kwargs['dt'] = self.dt
        super().__init__(self.dim, **kwargs)

        # Set the a0
        self.a0 = a0
        self.bool_has_defined_length_scale = True

        if not bool_is_for_properties_calculation:
            self.conf_apply_distortion([[final_strain,0],[0,final_strain]], update_q_and_a_vectors=True)
            self.a0 = self.a0 * (1+final_strain)
        
        
    def calc_proto_amplitudes_conserved(self):
        """Calculate the proto-amplitudes for the system.
        
        This method finds the optimal amplitude values (A, B) that minimize the free energy
        of the phase field crystal system with conserved dynamics.
        
        Returns
        -------
        tuple
            A tuple containing:
            - A : float
                The first proto-amplitude.
            - B : float
                The second proto-amplitude.
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
        """Calculate the equations for the proto-amplitudes for conserved dynamics.
        
        Parameters
        ----------
        vars : array_like
            The proto-amplitudes for the system [A, B].
            
        Returns
        -------
        list
            The equations for the proto-amplitudes that need to be solved.
            When both equations equal zero, the amplitudes are at equilibrium.
        """
        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v
        
        A, B = vars
        eq1 = A*(r + 9*A**2*v + 18*B**2*v + 2*t*psi0 + 3*v*psi0**2 + 4*B*(t + 3*v*psi0))
        eq2 = 9*B**3*v + 2*A**2*(t + 3*v*psi0) + B*(r + 18*A**2*v + 2*t*psi0 + 3*v*psi0**2)

        return [eq1, eq2]

    def calc_free_energy_from_proto_amplitudes(self, psi0, A, B):
        """Calculate the free energy of the system from the proto-amplitudes.
        
        Parameters
        ----------
        psi0 : float
            The average value of psi.
        A : float
            The first proto-amplitude.
        B : float
            The second proto-amplitude.
            
        Returns
        -------
        float
            The free energy of the system.
        """
        r = self.r
        t = self.t
        v = self.v
        return (108*A**4*v + 108*B**4*v + psi0**2*(24 + 6*r + 4*t*psi0 + 3*v*psi0**2) + 24*B**2*(r + psi0*(2*t + 3*v*psi0)) + 24*A**2*(r + 18*B**2*v + 4*B*(t + 3*v*psi0) + psi0*(2*t + 3*v*psi0)))
    
    def calc_L_f(self):
        """Calculates the L operator in Fourier space.

        Parameters
        ----------
        None

        Returns
        -------
        np.ndarray
            The L operator in Fourier space.
        """
        k2 = self.calc_k2()
        return (1 - k2)*(2 - k2)

    def calc_L_sum_f(self):
        """Calculate the sum of the L operators in Fourier space. Needed for stress calculation functions.

        Returns
        -------
        np.ndarray
            The sum of the L operators in Fourier space.
        """
        k2 = self.calc_k2()
        return 3 - 2*k2