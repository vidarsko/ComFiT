import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import matplotlib.pyplot as plt

class LandauSystem(cf.BaseSystem):
    
    def __init__(self,dim, r, **kwargs):
        self.r = r
        super().__init__(dim, **kwargs)


import scipy as sp

# Make linear operator
def calc_omega_f(self):
    return -self.calc_k2() - self.r

# Add method to class
LandauSystem.calc_omega_f = calc_omega_f

# Make the non-linear operator
def calc_nonlinear_evolution_function_f(self, field, t):
    return -sp.fft.fftn(field**3)

# Add method to class
LandauSystem.calc_nonlinear_evolution_function_f = calc_nonlinear_evolution_function_f

import numpy as np
def evolve(self, number_steps):
    omega_f = calc_omega_f(self)

    integrating_factors_f, solver = self.calc_integrating_factors_f_and_solver(omega_f, 'ETD2RK')

    for n in range(number_steps):
        self.psi, self.psi_f = solver(integrating_factors_f, 
                                    self.calc_nonlinear_evolution_function_f, 
                                    self.psi, self.psi_f)
        self.psi = np.real(self.psi)

# Add evolve method to class
LandauSystem.evolve = evolve

# Initiating a system with negative r value
# Decreasing the resolution a bit for faster computation
ls = LandauSystem(3, -0.5, xRes=31, yRes=31, zRes=31)

# Setting an initial condition with both positive and negative values
ls.psi = np.random.rand(ls.xRes,ls.yRes,ls.zRes)-0.5
ls.psi_f = sp.fft.fftn(ls.psi)
ls.evolve(300)
ls.plot_field(ls.psi,number_of_layers=5)
plt.show()