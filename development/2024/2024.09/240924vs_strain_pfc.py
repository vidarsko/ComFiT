import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


pfc = cf.PhaseFieldCrystal2DTriangular(5,2)
pfc.plot_lib = 'plotly'

pfc.conf_PFC_from_amplitudes()

distortion = [[0,0.6],
              [0.0,0]]

pfc.conf_apply_distortion(distortion)
# fig = pfc.plot_field(pfc.psi)
# fig.show()

ax, fig = pfc.plot_PFC()

alpha = pfc.calc_dislocation_density()
fig, ax1 = pfc.plot_field(alpha[0]**2+alpha[1]**2)
plt.show()
# plt.show()


# fig = pfc.plot_field(alpha[0]**2+alpha[1]**2)
# fig.show()
# for n in range(10):
#     pfc.evolve_PFC(10)
#     fig = pfc.plot_field(pfc.psi)
#     pfc.plot_save(fig, n)

# cf.tool_make_animation_gif(n)
