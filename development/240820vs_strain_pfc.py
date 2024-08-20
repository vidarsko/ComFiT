import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


pfc = cf.PhaseFieldCrystal2DTriangular(5,2)
pfc.plot_lib = 'plotly'

pfc.conf_PFC_from_amplitudes()
# fig = pfc.plot_field(pfc.psi)

for n in range(10):
    pfc.evolve_PFC(10)
    fig = pfc.plot_field(pfc.psi)
    cf.tool_save_plot_plotly(n,fig)

cf.tool_make_animation_gif(n)
