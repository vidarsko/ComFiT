from comfit.models.phase_field_crystal import PhaseFieldCrystal

import numpy as np
from tqdm import tqdm
from scipy.optimize import fsolve
import scipy as sp
import matplotlib.pyplot as plt
from pprint import pprint

class PhaseFieldCrystal2DTriangular(PhaseFieldCrystal):
    def __init__(self, nx, ny, **kwargs):
        """
        Nothing here yet
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

        self.A = self.calc_proto_amplitudes_conserved()
        self.eta0 = np.array([self.A, self.A, self.A])

        # Set the elastic constants
        self.el_lambda = 3 * self.A ** 2
        self.el_mu = 3 * self.A ** 2
        self.el_gamma = 0


        # Initialize the BaseSystem
        super().__init__(self.dim, xRes=self.xRes, yRes=self.yRes,
                         dx=self.dx, dy=self.dy, dt=self.dt)

        # Set the a0
        self.a0 = a0
        self.defined_length_scale = True

        

    def calc_proto_amplitudes_conserved(self):
        psi0 = self.psi0
        r = self.r
        t = self.t
        v = self.v
        
        #Testing which of the three amplitudes give the lowest free energy
        A = 0
        A_tmp = 1 / (15 * v)*(-t -3 * v * psi0 - np.sqrt(t**2 - 15 * v * r - 24 * t * v * psi0 - 36 * v ** 2 * psi0 ** 2)) 
        if self.calc_free_energy_from_proto_amplitudes(psi0, A_tmp) < self.calc_free_energy_from_proto_amplitudes(psi0, A):
            A = A_tmp
        A_tmp = 1 / (15 * v)*(-t -3 * v * psi0 + np.sqrt(t**2 - 15 * v * r - 24 * t * v * psi0 - 36 * v ** 2 * psi0 ** 2)) 
        if self.calc_free_energy_from_proto_amplitudes(psi0, A_tmp) < self.calc_free_energy_from_proto_amplitudes(psi0, A):
            A = A_tmp
        return A

    def calc_free_energy_from_proto_amplitudes(self, psi0, A):
        r = self.r
        t = self.t
        v = self.v

        return np.pi**2/(3*np.sqrt(3))*(270*A**4*v + 48*A**3*( t + 3*v*psi0) + psi0**2*(6 + 6*r + 4*t*psi0 + 3*v*psi0**2) + 36*A**2*( r + psi0*(2*t + 3*v*psi0)))

    def calc_omega_f(self):
        k2 = self.calc_k2()
        return - k2 * (self.r + (1 - k2)**2)
    
    def calc_stress_tensor_microscopic(self):
        """
        Calculates the microscopic stress of the phase-field crystal.
        """
        stress = np.zeros((3,self.xRes,self.yRes))
        
        Lpsi = np.real(sp.fft.ifftn((1-self.calc_k2())*self.psi_f))
        stress[0] = -2*Lpsi*np.real(sp.fft.ifftn(self.dif[0]*self.dif[0]*self.psi_f))
        stress[1] = -2*Lpsi*np.real(sp.fft.ifftn(self.dif[0]*self.dif[1]*self.psi_f))
        stress[2] = -2*Lpsi*np.real(sp.fft.ifftn(self.dif[1]*self.dif[1]*self.psi_f))

        return stress

    
    def calc_stress_divergence_f(self, field_f = None):
        
        if field_f is None:
            field_f = self.psi_f
        
        k2 = self.calc_k2()

        return np.array([
            -2*self.calc_Gaussian_filter_f()*sp.fft.fftn(sum([
                sp.fft.ifftn((1-k2)*self.dif[i]*field_f)*sp.fft.ifftn(self.dif[i]*self.dif[j]*field_f) for i in range(self.dim)
                ]) + sp.fft.ifftn((1-k2)*field_f)*sp.fft.ifftn(self.dif[j]*(-k2)*field_f)) for j in range(self.dim)]
                )