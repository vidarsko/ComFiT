import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


size = 20
pfc = cf.PhaseFieldCrystal2DTriangular(size,round(size/np.sqrt(3)),type_of_evolution='', t=1/2, r=0.02)       
print(pfc.psi0)
print(pfc.A)                            
# pfc.conf_PFC_from_amplitudes()
# pfc.evolve_PFC


# distortion = [[0.3,0.0],
#               [0.0,0]]

# pfc.conf_apply_distortion(distortion)
# # fig = pfc.plot_field(pfc.psi)

# # fig.show()

# noise_strength = 0.01
# for n in range(100):
#     pfc.psi = pfc.psi + noise_strength*np.random.randn(*pfc.psi.shape)
#     pfc.psi_f = sp.fft.fftn(pfc.psi)
#     pfc.evolve_PFC_unconserved(10)
#     fig, axs = plt.subplots(1,2)
#     pfc.plot_PFC(ax=axs[0],title='PFC')

#     alpha = pfc.calc_dislocation_density()
#     pfc.plot_field(alpha[0]**2+alpha[1]**2, ax=axs[1], title='Dislocation density')
#     cf.tool_save_plot_matplotlib(n)
# cf.tool_make_animation_gif(n)

