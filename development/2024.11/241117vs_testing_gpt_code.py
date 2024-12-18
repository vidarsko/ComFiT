import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


# Initialize a Bose-Einstein Condensate (BEC) system
# bec = cf.BoseEinsteinCondensate(dim=2, xRes=256, yRes=256, xmin=-20, xmax=20, ymin=-20, ymax=20, dt=0.01)
# bec = cf.BoseEinsteinCondensate(dim=2)
# print(bec.dx)
# print(bec.dy)
# Configure the initial condition with a vortex-antivortex pair
# vortex_charge = 1  # Positive vortex
# antivortex_charge = -1  # Negative vortex
# vortex_position = [5, 0]  # Position of vortex
# antivortex_position = [-5, 0]  # Position of antivortex
# bec.conf_insert_vortex()
# bec.conf_insert_vortex(charge=vortex_charge, position=vortex_position)
# bec.conf_insert_vortex(charge=antivortex_charge, position=antivortex_position)

# Relax the system to dissipate any unwanted artifacts
# bec.evolve_relax(0)

# Simulate vortex-antivortex annihilation
# bec.evolve_dGPE(number_of_steps=0)

# Plot the resulting field
# fig = bec.plot_complex_field(bec.psi)
# fig.show()


# bs = cf.BaseSystem(2, xmin=0, xmax=20, xlim=[-2,22], dx = 3, xRes=256, ylim = [0, 20], ymax = 20, dy = 3, yRes = 256)


# bec = cf.BoseEinsteinCondensate(2,xRes=31,yRes=31)
# print(bec.xmax)
# bec.conf_insert_vortex_dipole(
#     dipole_vector=[bec.xmax/3,0],
#     dipole_position=bec.rmid)
# bec.evolve_relax(100)
# Dnodes = bec.calc_vortex_nodes()

# import comfit as cf
# import numpy as np
# import scipy as sp

class CahnHilliardSystem(cf.BaseSystem):
    def __init__(self, dim, gamma, **kwargs):
        self.gamma = gamma
        super().__init__(dim, **kwargs)
    
    def calc_nonlinear_evolution_function_f(self, field, t):
        # Nonlinear term: \nabla^2(\psi^3 - \psi)
        laplacian_f = -self.calc_k2()  # \nabla^2 operator in Fourier space
        psi3_minus_psi_f = sp.fft.fftn(field**3 - field)  # Transform (\psi^3 - \psi) to Fourier space
        return laplacian_f * psi3_minus_psi_f
    
    def evolve(self, number_steps):
        # Conserved evolution
        omega_f = -self.gamma * self.calc_k2()**2  # Linear operator: -γ \nabla^4
        integrating_factors_f, solver = self.calc_integrating_factors_f_and_solver(omega_f, method='ETD4RK')
        for n in range(number_steps):
            self.psi, self.psi_f = solver(integrating_factors_f, 
                                          self.calc_nonlinear_evolution_function_f, 
                                          self.psi, self.psi_f)
            self.psi = np.real(self.psi)  # Keep \psi real since it's real-valued

# Example simulation
gamma = 1.0
ch = CahnHilliardSystem(dim=2, gamma=gamma, xRes=128, yRes=128, xmin=0, xmax=10, ymin=0, ymax=10, dt=0.01)
ch.psi = np.random.rand(ch.xRes, ch.yRes) - 0.5  # Random initial condition
ch.psi_f = sp.fft.fftn(ch.psi)  # Fourier transform of initial condition

# Evolve system
ch.evolve(200)

# Plot result
ch.plot_field(ch.psi).show()
