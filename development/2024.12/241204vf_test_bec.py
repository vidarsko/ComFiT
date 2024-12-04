import comfit as cf
import plotly
import numpy as np


bec = cf.BoseEinsteinCondensate(2)
bec.conf_insert_vortex_dipole()

bec.evolve_relax(100)
fig1 = bec.plot_complex_field(bec.psi)
fig2 = bec.plot_field(abs(bec.psi))
fig = cf.tool_make_subplots(1, 2, fig1, fig2)
fig.show()
