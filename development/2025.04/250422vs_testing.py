import comfit as cf

pfc = cf.PhaseFieldCrystal2DSquare(30,30)
eta = pfc.calc_amplitudes_with_dislocation_dipole(dislocation_type=2)
# eta = pfc.calc_amplitudes_with_dislocation(dislocation_type=1)
eta = pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(100)

fig, ax = pfc.plot_field(pfc.psi)
pfc.show(fig)

import matplotlib.pyplot as plt
import numpy as np

psi0 = pfc.ifft(pfc.psi_f*pfc.calc_Gaussian_filter_f()).real

dxpsi0 = pfc.ifft(pfc.dif[0]*pfc.fft(psi0)).real
dypsi0 = pfc.ifft(pfc.dif[1]*pfc.fft(psi0)).real

abs_val = dxpsi0**2 + dypsi0**2

fig, ax = pfc.plot_field(abs_val)
pfc.show(fig)

# # Plot a cross-section of psi at the middle y value
# plt.figure()
# plt.plot(pfc.x, np.gradient(psi0[:, pfc.ymidi],pfc.dx))
# plt.xlabel('x')
# plt.ylabel('psi')
# plt.title('Cross-section of psi at y = middle')
# plt.grid(True)
# plt.show()