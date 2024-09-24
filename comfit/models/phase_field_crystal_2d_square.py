from comfit.models.phase_field_crystal import PhaseFieldCrystal

import numpy as np
from tqdm import tqdm
from scipy.optimize import fsolve
import scipy as sp
import matplotlib.pyplot as plt
from pprint import pprint

class PhaseFieldCrystal2DSquare(PhaseFieldCrystal):
    def __init__(self, nx, ny, **kwargs):
        """Initializes a phase field crystal system in 2D with a square crystal structure.

        Args:
            nx: The number of unit cells in the x direction.
            ny: The number of unit cells in the y direction.

        Returns:
            The system object representing the PhaseFieldCrystal2DSquare simulation.
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
        self.q = np.array([[1, 0], [0, 1], [1, -1], [1, 1]])

        # Set the number of reciprocal modes
        self.number_of_reciprocal_lattice_modes = 4
        self.number_of_primary_reciprocal_lattice_modes = 2

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]

        self.A, self.B = self.calc_proto_amplitudes_conserved()
        self.eta0 = np.array([self.A, self.A, self.B, self.B])

        # Set the elastic constants
        self.el_lambda = 16 * self.B ** 2
        self.el_mu = 16 * self.B ** 2
        self.el_gamma = 0


        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes,
                         dx=self.dx, dy=self.dy, dt=self.dt)

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

        A = 0
        B = 0
        free_energy = self.calc_free_energy_from_proto_amplitudes(self.psi0, A, B)

        for A0 in np.linspace(0.01,1,10):
            for B0 in [A0/2]:
                [A_tmp, B_tmp] = fsolve(self.calc_proto_amplitude_equations_conserved, [A0, B0])
                free_energy_tmp = self.calc_free_energy_from_proto_amplitudes(self.psi0, A_tmp, B_tmp)

                if free_energy_tmp <= free_energy:
                    A = A_tmp
                    B = B_tmp
                    free_energy = free_energy_tmp

        return A, B

    def calc_proto_amplitude_equations_conserved(self,vars):
        """Calculates the equations for the proto-amplitudes for the system in the case of conserved dynamics

        Args:
            vars: The proto-amplitudes for the system.

        Returns:
            The equations for the proto-amplitudes for the system.
        """
        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v
        
        A, B = vars
        eq1 = A*(r + 9*A**2*v + 18*B**2*v + 2*t*psi0 + 3*v*psi0**2 + 4*B*(t + 3*v*psi0))
        eq2 = 9*B**3*v + 2*A**2*(t + 3*v*psi0) + B*(r + 18*A**2*v + 2*t*psi0 + 3*v*psi0**2)

        return [eq1, eq2]

    def calc_free_energy_from_proto_amplitudes(self,psi0,A,B):
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
        return (108*A**4*v + 108*B**4*v + psi0**2*(24 + 6*r + 4*t*psi0 + 3*v*psi0**2) + 24*B**2*(r + psi0*(2*t + 3*v*psi0)) + 24*A**2*(r + 18*B**2*v + 4*B*(t + 3*v*psi0) + psi0*(2*t*+ 3*v*psi0)))

    def calc_chemical_potential_linear_part_f(self):
        """Calcuates the linear part of the chemical potential of the system.

        Args:
            None
        Returns:
            The linear part of the chemical potential of the system.
        """
        k2 = self.calc_k2()
        return (self.r + (1 - k2)**2*(2-k2)**2)

    def calc_omega_f(self):
        """Calculates the linear time evolution operator in Fourier space.

        Args:
            None

        Returns:
            The linear time evolution operator in Fourier space.
        """
        k2 = self.calc_k2()
        return - k2 * (self.r + (1 - k2)**2*(2-k2)**2)

    def calc_stress_tensor_microscopic(self):
        """Calculates the microscopic stress of the phase-field crystal.

        Args:
            None

        Returns:
            The microscopic stress of the phase-field crystal.
        """
        stress = np.zeros((3,self.xRes,self.yRes))
        
        k2 = self.calc_k2()
        L1L2psi = np.real(sp.fft.ifftn((1-k2)*(2-k2)*self.psi_f))
        L1plusL2_f = 2 - 3*k2
        stress[0] = -2*L1L2psi*np.real(sp.fft.ifftn(L1plusL2_f*self.dif[0]*self.dif[0]*self.psi_f))
        stress[1] = -2*L1L2psi*np.real(sp.fft.ifftn(L1plusL2_f*self.dif[0]*self.dif[1]*self.psi_f))
        stress[2] = -2*L1L2psi*np.real(sp.fft.ifftn(L1plusL2_f*self.dif[1]*self.dif[1]*self.psi_f))

        return stress

    def calc_stress_divergence_f(self, field_f = None):
        """Calculates the divergence of the stress tensor in Fourier space.
        
        Args:
            field_f: The field in Fourier space.
        
        Returns:
            The divergence of the stress tensor in Fourier space.
        """
        
        if field_f is None:
            field_f = self.psi_f
        
        k2 = self.calc_k2()

        return np.array([
            -2*self.calc_Gaussian_filter_f()*sp.fft.fftn(
                sum([
                sp.fft.ifftn((1-k2)*(2-k2)*self.dif[i]*field_f)*sp.fft.ifftn((3-2*k2)*self.dif[i]*self.dif[j]*field_f) 
                for i in range(self.dim)
                ]) 
                +sp.fft.ifftn((1-k2)*(2-k2)*field_f)*sp.fft.ifftn((3-2*k2)*self.dif[j]*(-k2)*field_f)) 
                for j in range(self.dim)]
                )
            #TODO: Needs double checking (Vidar 27.02.24)