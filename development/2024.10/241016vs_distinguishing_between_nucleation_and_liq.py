import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp
import copy


# Parameters
nx = 20
B_x = 0.98
psi0 = 0
t = 1 / 2 / B_x
v = 1 / 3 / B_x
time_limit = 100000


Delta_B = -0.02/B_x

pfc = cf.PhaseFieldCrystal2DTriangular(nx,round(nx/np.sqrt(3)), t=t, r=Delta_B, v=v, psi0=psi0)
print(pfc.A_proto)
pfc.conf_PFC_from_amplitudes()

strain=-0.3

distortion = np.eye(2)*strain
pfc.conf_apply_distortion(distortion, update_q_and_a_vectors=True)

pfc.plot_lib = 'matplotlib'

fig, axs = plt.subplots(1,2)


number_of_points = 100
points = np.arange(number_of_points)


eta_max = np.nan * np.zeros(number_of_points)
eta_min = np.nan * np.zeros(number_of_points)

noise_strength = 0.01

for n in points:
    
    pfc.psi = pfc.psi + noise_strength*np.random.randn(*pfc.psi.shape)
    pfc.psi_f = sp.fft.fftn(pfc.psi)    
    pfc.evolve_PFC(100)
    eta = pfc.calc_demodulate_PFC()
    eta_max[n] = np.max(abs(eta))
    eta_min[n] = np.min(abs(eta))


    # Clear the figure but keep the subplots
    for ax in axs:
        ax.clear()

    pfc.plot_field(pfc.psi, ax=axs[0], title='PFC', colorbar=False)

    # axs[1].plot(points, eta_max, label='eta_max')
    # axs[1].plot(points, eta_min, label='eta_min')
    axs[1].plot((eta_max-eta_min)/eta_max, label='relative difference')
    axs[1].set_xlabel('Points')
    axs[1].legend()
    plt.draw()
    plt.pause(0.01)
    
