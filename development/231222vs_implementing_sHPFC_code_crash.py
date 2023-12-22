import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

Delta_t = 2
dt = 0.1

pfc = cf.PhaseFieldCrystal2DTriangular(6,3, dt=dt)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)
# pfc.conf_PFC_from_amplitudes()


for n in range(40):
    plt.clf()
    pfc.evolve_PFC_hydrodynamic(round(Delta_t/dt), rho0=2**-2, gamma_S=2**-2)
    f = np.real(np.fft.ifftn(pfc.calc_stress_divergence_f(pfc.psi_f[0]), axes = (range(-pfc.dim,0))))
    fnorm = np.sqrt(f[0]**2 + f[1]**2)
    #pfc.plot_field(fnorm,colormap='winter')
    #pfc.plot_field(f[0],colormap='winter')

    #vnorm = np.sqrt(pfc.psi[1]**2 + pfc.psi[2]**2)
    # pfc.plot_field(vnorm, colormap='winter')
    #pfc.plot_vector_field([pfc.psi[1],pfc.psi[2]])

    # pfc.plot_angle_field(np.arctan2(pfc.psi[2],pfc.psi[1]),colorbar=False)
    pfc.plot_field(pfc.psi[0],colormap='winter')
    # pfc.plot_field(pfc.psi[1],colormap='winter')
    plt.draw()
    plt.pause(0.01)

plt.savefig(f'dt={dt}.png')