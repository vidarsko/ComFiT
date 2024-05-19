import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import plotly.graph_objects as go

bec = cf.BoseEinsteinCondensate(2,xRes=51,yRes=101)
bec.a0=3
# bec = cf.BoseEinsteinCondensate(3,xRes=31,yRes=31,zRes=31)
bec.plot_lib = 'plotly'
# bec.conf_insert_vortex_dipole()
# bec.conf_insert_vortex_ring()
# bec.evolve_relax(100)
bec.psi = bec.calc_Gaussian(width=10)*np.exp(1j*bec.x/3)
# fig = bec.plot_field(abs(bec.psi),number_of_layers=4,alpha=0.1)
fig = bec.plot_complex_field(bec.psi)
# fig = bec.plot_field(abs(bec.psi))
# bec.show()
fig.show()
# plt.show()