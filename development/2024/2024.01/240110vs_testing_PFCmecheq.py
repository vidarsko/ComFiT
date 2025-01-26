import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np


# 240110vs_testing_PFCmecheq.py

pfc = cf.PhaseFieldCrystal2DTriangular(21,12)
# pfc = cf.PhaseFieldCrystal2DSquare(21,21)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(100)

fig,axs = plt.subplots(2,2)

pfc.plot_field(pfc.psi,ax=axs[0,0],colorbar=False)

# u = pfc.calc_displacement_field_to_equilibrium()
u = np.zeros((2,pfc.xRes,pfc.yRes))
u[0] = np.sin(pfc.y/pfc.ymax*2*np.pi)

g_f = pfc.calc_stress_divergence_f()
g = np.real(np.fft.ifftn(g_f,  axes = (range (-2, 0) ) ))

pfc.conf_advect_PFC(u) #seems like this is ok. 
# print(pfc.psi)

pfc.plot_vector_field(u,ax=axs[0,1])


pfc.plot_vector_field(g,ax=axs[1,1])

pfc.plot_field(pfc.psi,ax=axs[1,0],colorbar=False)
plt.show()

