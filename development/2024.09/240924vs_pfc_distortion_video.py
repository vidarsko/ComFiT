import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp

pfc = cf.PhaseFieldCrystal2DTriangular(20,round(20/np.sqrt(3)))
print(pfc.psi0)
print(pfc.A)
pfc.conf_PFC_from_amplitudes()

distortion = [[0.0,0.15],
              [0.15,0]]

pfc.conf_apply_distortion(distortion)
# fig = pfc.plot_field(pfc.psi)

# fig.show()

noise_strength = 0.01
for n in range(100):
    pfc.psi = pfc.psi + noise_strength*np.random.randn(*pfc.psi.shape)
    pfc.psi_f = sp.fft.fftn(pfc.psi)
    pfc.evolve_PFC(100)
    fig, axs = plt.subplots(1,2)
    pfc.plot_PFC(ax=axs[0],title='PFC')

    alpha = pfc.calc_dislocation_density()
    pfc.plot_field(alpha[0]**2+alpha[1]**2, ax=axs[1], title='Dislocation density')
    cf.tool_save_plot_matplotlib(n)
cf.tool_make_animation_gif(n)

