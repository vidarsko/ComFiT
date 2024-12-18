import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


bec = cf.BoseEinsteinCondensate(2)
bec.conf_insert_vortex_dipole()

bec.evolve_relax(100)
fig1 = bec.plot_complex_field(bec.psi)
fig2 = bec.plot_field(abs(bec.psi))
fig = cf.tool_make_subplots(1, 2, fig1, fig2)
fig.show()

cf.tool_make_animation_gif
