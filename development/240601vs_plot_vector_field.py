import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import time

import plotly.graph_objects as go

nx = 12
pfc = cf.PhaseFieldCrystal2DTriangular(nx,np.round(nx/np.sqrt(3)))
pfc.plot_lib = 'plotly'

pfc.conf_PFC_from_amplitudes()
gradpsi = np.array([sp.fft.ifftn(pfc.dif[0]*pfc.psi_f), sp.fft.ifftn(pfc.dif[1]*pfc.psi_f)])
# gradpsi = np.array([1*(pfc.x+pfc.y), 1*(pfc.x-pfc.y)])

# fig = go.Figure(data=[go.Scatter(x=pfc.x, y=pfc.y, mode='markers', marker=dict(size=5, color=pfc.psi, colorscale='Viridis', showscale=True))])

# fig.show()

start = time.time()
fig = pfc.plot_vector_field(gradpsi, spacing=1)
# plt.show()
fig.show()
print(time.time()-start)