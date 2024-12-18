import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


pfc = cf.PhaseFieldCrystal2DSquare(20,20)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)

pfc.evolve_PFC(100)

nodes = pfc.calc_dislocation_nodes()
print(nodes)

fig = pfc.plot_PFC()
pfc.plot_nodes(nodes,fig=fig)

# fig = pfc.plot_nodes(nodes)

fig.show()