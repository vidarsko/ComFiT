import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

bec = cf.BoseEinsteinCondensate(2, xRes=51, yRes=51, gamma=0)
# bec = cf.BoseEinsteinCondensate(2, xRes=25, yRes=25, gamma=0)

# bec.conf_insert_vortex_dipole(dipole_vector=[14,0])
bec.conf_insert_vortex_dipole()
bec.evolve_relax(200)
roll = 30

# dens = bec.calc_vortex_density()
# bec.plot_field(abs(dens),colormap='bluewhitered',vlim_symmetric=True)
# plt.show()

for n in range(30):
    bec.evolve_dGPE(100)
    plt.clf()
    _, ax = bec.plot_complex_field(bec.psi)
    psi_old = bec.psi.copy()
    bec.evolve_dGPE(10)
    dt_psi = (bec.psi - psi_old)/(10*bec.dt)
    Dnodes = bec.calc_vortex_nodes(dt_psi=dt_psi)
    bec.plot_nodes(Dnodes, ax=ax)
    plt.pause(0.01)
plt.show()