import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


bec = cf.BoseEinsteinCondensate(3, xRes=31, yRes=31, zRes=31)
bec.conf_insert_vortex_ring()
bec.evolve_relax(300)

# fig = bec.plot_field(abs(bec.psi))
fig = bec.plot_complex_field(bec.psi, phase_blob_threshold=0.9, interpolation_method='nearest')
# fig = bec.plot_complex_field_in_plane(bec.psi)

fig.show()
