from comfit.models.phase_field_crystal import PhaseFieldCrystal

import numpy as np
from tqdm import tqdm
from scipy.optimize import fsolve
import scipy as sp
import matplotlib.pyplot as plt
from pprint import pprint

class PhaseFieldCrystal3DSimpleCubic(PhaseFieldCrystal):
    def __init__(self, nx, ny, nz, **kwargs):
        """Initializes a phase field crystal system in 3D with a simple cubic crystal structure.

        Args:
            nx: The number of unit cells in the x direction.
            ny: The number of unit cells in the y direction.
            nz: The number of unit cells in the z direction.

        Returns:
            The system object representing the PhaseFieldCrystal3DSimpleCubic simulation.
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
                           [1, 1, 1]])

        # Set the number of reciprocal modes
        self.number_of_reciprocal_lattice_modes = 13
        self.number_of_primary_reciprocal_lattice_modes = 3

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]
        self.dz = a0 / self.micro_resolution[2]

        self.A, self.B, self.C = self.calc_proto_amplitudes_conserved()
        self.eta0 = np.array([self.A, self.A, self.A,
                     self.B, self.B, self.B, self.B, self.B, self.B,
                     self.C, self.C, self.C, self.C])

        # Set the elastic constants
        self.el_lambda = 16 * self.B ** 2 + 128 * self.C ** 2
        self.el_mu = 16 * self.B ** 2 + 128 * self.C ** 2
        self.el_gamma = 32*self.A**2 - 16*self.B**2 - 256*self.C**2

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes, zRes=self.zRes,
                         dx=self.dx, dy=self.dy, dz=self.dz, dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

    def calc_free_energy_from_proto_amplitudes(self, psi0, A, B, C):
        """Calculates the free energy of the system from the proto-amplitudes.

        Args:
            psi0: The value of psi0.
            A: The value of A.
            B: The value of B.
            C: The value of C.

        Returns:
            The free energy of the system.
        """

        r = self.r
        t = self.t
        v = self.v
        return 2*np.pi**3/3*(48*C**2*r + 270*A**4*v + 1620*B**4*v + 576*A**3*C*v + 648*C**4*v + 96*C**2*t*psi0 + 6*(36 + r + 24*C**2*v)*psi0**2 + 4*t*psi0**3 + 3*v*psi0**4 + 192*B**3*(t + 3*v*psi0) + 576*A*B*C*(t + 3*v*(3*B + psi0)) + 36*A**2*( r + 96*B**2*v + 36*C**2*v + 2*t*psi0 + 3*v*psi0**2 + 8*B*( t + 3*v*psi0)) + 72*B**2*(r + 54*C**2*v + psi0*(2*t + 3*v*psi0)) )

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

        A, B, C = vars

        eq1 = 15*A**3*v + 24*A**2*C*v + 8*B*C*(t + 3*v*(3*B + psi0)) + A*( r + 96*B**2*v + 36*C**2*v + 2*t*psi0 + 3*v*psi0**2 + 8*B*(t + 3*v*psi0))
        eq2 = 45*B**3*v + 4*B**2*(t + 3*v*psi0) + 2*A*(A + 2*C)*(t + 3*v*psi0) + B*( r + 48*A**2*v + 72*A*C*v + 54*C**2*v + 2*t*psi0 + 3*v*psi0**2)
        eq3 = 27*C**3*v + C*( r + 27*A**2*v + 81*B**2*v + 2*t*psi0 + 3*v*psi0**2) + 6*A*(A**2*v + 9*B**2*v + B*(t + 3*v*psi0))

        return [eq1, eq2, eq3]

    def calc_L_f(self):
        """Calculates the L operator in Fourier space.

        Args:
            None
        Returns:
            The L operator in Fourier space.
        """
        k2 = self.calc_k2()
        return (1 - k2)*(2 - k2)*(3 - k2)

    def calc_stress_tensor_microscopic(self):
        """Calculates the microscopic stress of the phase-field crystal.

        Args:
            None

        Returns:
            The microscopic stress of the phase-field crystal.
        """
        stress = np.zeros((6,self.xRes,self.yRes,self.zRes))
        
        k2 = self.calc_k2()
        L1L2L3psi = np.real(sp.fft.ifftn((1-k2)*(2-k2)*(3-k2)*self.psi_f))
        L2L3_p_L1L3_p_L1L2_f = 11 - 12*k2 + 3*k2**2
        stress[0] = -2*L1L2L3psi*np.real(sp.fft.ifftn(L2L3_p_L1L3_p_L1L2_f*self.dif[0]*self.dif[0]*self.psi_f))
        stress[1] = -2*L1L2L3psi*np.real(sp.fft.ifftn(L2L3_p_L1L3_p_L1L2_f*self.dif[0]*self.dif[1]*self.psi_f))
        stress[2] = -2*L1L2L3psi*np.real(sp.fft.ifftn(L2L3_p_L1L3_p_L1L2_f*self.dif[0]*self.dif[2]*self.psi_f))
        stress[3] = -2*L1L2L3psi*np.real(sp.fft.ifftn(L2L3_p_L1L3_p_L1L2_f*self.dif[1]*self.dif[1]*self.psi_f))
        stress[4] = -2*L1L2L3psi*np.real(sp.fft.ifftn(L2L3_p_L1L3_p_L1L2_f*self.dif[1]*self.dif[2]*self.psi_f))
        stress[5] = -2*L1L2L3psi*np.real(sp.fft.ifftn(L2L3_p_L1L3_p_L1L2_f*self.dif[2]*self.dif[2]*self.psi_f))

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
            -2*self.calc_Gaussian_filter_f()*sp.fft.fftn(sum([
                sp.fft.ifftn((1-k2)*(2-k2)*(3-k2)*self.dif[i]*field_f)*sp.fft.ifftn((11-12*k2+3*k2**2)*self.dif[i]*self.dif[j]*field_f) for i in range(self.dim)
                ]) + sp.fft.ifftn((1-k2)*(2-k2)*(3-k2)*field_f)*sp.fft.ifftn((11-12*k2+3*k2**2)*self.dif[j]*(-k2)*field_f)) for j in range(self.dim)]
                )