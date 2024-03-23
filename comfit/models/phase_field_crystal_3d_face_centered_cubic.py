from comfit.models.phase_field_crystal import PhaseFieldCrystal

import numpy as np
from tqdm import tqdm
from scipy.optimize import fsolve
import scipy as sp
import matplotlib.pyplot as plt
from pprint import pprint

class PhaseFieldCrystal3DFaceCenteredCubic(PhaseFieldCrystal):
    def __init__(self, nx, ny, nz, **kwargs):
        """
        Nothing here yet
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
                                    [-1, 1, 0]])

        self.q = np.array([[-1, 1, 1],
                           [1, -1, 1],
                           [1, 1, -1],
                           [1, 1, 1],
                           [2, 0, 0],
                           [0, 2, 0],
                           [0, 0, 2]]) / np.sqrt(3)

        # Set the number of reciprocal modes
        self.number_of_reciprocal_lattice_modes = 7
        self.number_of_primary_reciprocal_lattice_modes = 4

        # Set the grid
        self.dx = a0 / self.micro_resolution[0]
        self.dy = a0 / self.micro_resolution[1]
        self.dz = a0 / self.micro_resolution[2]

        self.A, self.B = self.calc_proto_amplitudes_conserved()
        self.eta0 = np.array([self.A, self.A, self.A, self.A, self.B, self.B, self.B])

        # Set the elastic constants
        self.el_lambda = 32/81 * self.A ** 2
        self.el_mu = 32/81 * self.A ** 2
        self.el_gamma = 32/81 * (2*self.B**2 - self.A**2)

        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes, zRes=self.zRes,
                         dx=self.dx, dy=self.dy, dz=self.dz, dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

        
    def calc_proto_amplitudes_conserved(self):
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

        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v

        A, B = vars
        eq1 = A*( r + 27*A**2*v + 36*B**2*v + 2*t*psi0 + 3*v*psi0**2 + 6*B*(t + 3*v*psi0))
        eq2 = 15*B**3*v + 4*A**2*(t + 3*v*psi0) + B*( r + 48*A**2*v + 2*t*psi0 + 3*v*psi0**2)

        return [eq1, eq2]

    def calc_free_energy_from_proto_amplitudes(self, psi0, A, B):
        r = self.r
        t = self.t
        v = self.v
        return 2*np.pi**3/np.sqrt(3)*(1944*A**4*v + 810*B**4*v + psi0**2*(32 + 18*r + 12*t*psi0 + 9*v*psi0**2) + 108*B**2*(r + psi0*(2*t + 3*v*psi0)) + 144*A**2*( r + 36*B**2*v + 6*B*(t + 3*v*psi0) + psi0*(2*t + 3*v*psi0)))

    def calc_omega_f(self):
        k2 = self.calc_k2()
        return - k2 * (self.r + (1 - k2)**2*(4/3-k2)**2)

    def calc_stress_tensor_microscopic(self):
        """
        Calculates the microscopic stress of the phase-field crystal.
        """
        stress = np.zeros((6,self.xRes,self.yRes,self.zRes))
        
        k2 = self.calc_k2()
        L1L43psi = np.real(sp.fft.ifftn((1-k2)*(4/3-k2)*self.psi_f))
        L1plusL43_f = 7/3 - 2*k2
        stress[0] = -2*L1L43psi*np.real(sp.fft.ifftn(L1plusL43_f*self.dif[0]*self.dif[0]*self.psi_f))
        stress[1] = -2*L1L43psi*np.real(sp.fft.ifftn(L1plusL43_f*self.dif[0]*self.dif[1]*self.psi_f))
        stress[2] = -2*L1L43psi*np.real(sp.fft.ifftn(L1plusL43_f*self.dif[0]*self.dif[2]*self.psi_f))
        stress[3] = -2*L1L43psi*np.real(sp.fft.ifftn(L1plusL43_f*self.dif[1]*self.dif[1]*self.psi_f))
        stress[4] = -2*L1L43psi*np.real(sp.fft.ifftn(L1plusL43_f*self.dif[1]*self.dif[2]*self.psi_f))
        stress[5] = -2*L1L43psi*np.real(sp.fft.ifftn(L1plusL43_f*self.dif[2]*self.dif[2]*self.psi_f))

        return stress

    def calc_stress_divergence_f(self, field_f = None):
        
        if field_f is None:
            field_f = self.psi_f
        
        k2 = self.calc_k2()

        return np.array([
            -2*self.calc_Gaussian_filter_f()*sp.fft.fftn(sum([
                sp.fft.ifftn((1-k2)*(4/3-k2)*self.dif[i]*field_f)*sp.fft.ifftn((7/3 - 2*k2)*self.dif[i]*self.dif[j]*field_f) for i in range(self.dim)
                ]) + sp.fft.ifftn((1-k2)*(4/3-k2)*field_f)*sp.fft.ifftn(7/3 - 2*k2)*self.dif[j]*(-k2)*field_f) for j in range(self.dim)])
                #TODO: Needs double checking (Vidar 27.02.24)