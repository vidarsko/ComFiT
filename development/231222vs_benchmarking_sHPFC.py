import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

pfc = cf.PhaseFieldCrystal2DTriangular(12,6)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(100)
f = np.real(np.fft.ifftn(pfc.calc_stress_divergence_f(pfc.psi_f), axes = (range(-pfc.dim,0))))
fnorm = np.sqrt(f[0]**2 + f[1]**2)/(pfc.el_mu/pfc.a0)
pfc.plot_field(fnorm,colormap='parula',cmap_symmetric=False)
# pfc.plot_field(pfc.psi,colormap='parula')
print("mu",pfc.el_mu)
print("a0",pfc.a0)
print(pfc.A)
plt.show()
