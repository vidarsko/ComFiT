import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

bec = cf.BoseEinsteinCondensate(2)
bec.conf_initial_condition_disordered()

# bec.V_ext = 0

# def V_ext(self):
#     psi_ideal = 1
#     return abs(self.psi-psi_ideal)**2

# bec.V_ext = V_ext.__get__(bec)

def calc_nonlinear_evolution_function_f(self, psi):
        psi2 = np.abs(psi) ** 2

        A=10
        psi_ideal = np.exp(1j*np.arctan2(self.y-self.ymid,self.x-self.xmid))
        radius = self.xmax/2
        mask = (self.x-self.xmid)**2 + (self.y-self.ymid)**2 > (radius)**2


        return sp.fft.fftn((1j + self.gamma) * (-self.V_ext() - psi2) * psi - mask*A*(psi-psi_ideal))

bec.calc_nonlinear_evolution_function_f = calc_nonlinear_evolution_function_f.__get__(bec)

for n in range(20):
    plt.clf()
    bec.evolve_relax(200)
    bec.plot_complex_field(bec.psi)
    # bec.plot_field(abs(bec.psi),colorbar=True)
    plt.draw()
    plt.pause(0.01)

plt.show()
