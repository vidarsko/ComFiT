import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp

nx = 20
pfc = cf.PhaseFieldCrystal2DTriangular(nx,round(nx/np.sqrt(3)))
# pfc = cf.PhaseFieldCrystal2DSquare(20,20)
pfc.plot_lib = 'matplotlib'
pfc.conf_PFC_from_amplitudes()

distortion = [[0.0,0.15],
              [0.15,0.0]]

pfc.conf_apply_distortion(distortion)
print(pfc.psi)
# pfc.evolve_PFC(100)

def plot():
    fig, axs = plt.subplots(1,2)
    pfc.plot_PFC(ax=axs[0],title='PFC')

    alpha = pfc.calc_dislocation_density()
    pfc.plot_field(np.sqrt(alpha[0]**2+alpha[1]**2), vlim=[0,0.0045], ax=axs[1], title='Dislocation density')

# plot()
# plt.show()

# fig = pfc.plot_field(pfc.psi)

# fig.show()

noise_strength = 0.01
max_alpha = 0
for n in range(100):
    pfc.psi = pfc.psi + noise_strength*np.random.randn(*pfc.psi.shape)
    pfc.psi_f = sp.fft.fftn(pfc.psi)
    pfc.evolve_PFC(100)
    alpha = pfc.calc_dislocation_density()
    max_alpha = max(max_alpha, np.sqrt(alpha[0]**2+alpha[1]**2).max())
    plot()
    cf.tool_save_plot_matplotlib(n)
# print(max_alpha)
cf.tool_make_animation_gif(n)

